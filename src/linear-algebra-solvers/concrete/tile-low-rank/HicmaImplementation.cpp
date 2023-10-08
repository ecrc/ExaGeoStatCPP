
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file HicmaImplementation.cpp
 * @brief Sets up the HiCMA descriptors needed for the tile low rank computations in ExaGeoStat.
 * @version 1.0.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-03-26
**/

#include <linear-algebra-solvers/concrete/hicma/tile-low-rank/HicmaImplementation.hpp>

using namespace std;

using namespace exageostat::linearAlgebra::tileLowRank;
using namespace exageostat::common;
using namespace exageostat::dataunits;
using namespace exageostat::kernels;
using namespace exageostat::helpers;
using namespace exageostat::hardware;
using namespace exageostat::configurations;

//// TODO: These variables are required to avoid undefined reference to HiCMA global variables.
int store_only_diagonal_tiles = 1;
int use_scratch = 1;
int global_check = 0;  //used to create dense matrix for accuracy check

template<typename T>
void HicmaImplementation<T>::SetModelingDescriptors(ExaGeoStatData<T> &aData, Configurations &aConfigurations) {

    int N = aConfigurations.GetProblemSize();
    int lts = aConfigurations.GetLowTileSize();
    int p_grid = aConfigurations.GetPGrid();
    int q_grid = aConfigurations.GetQGrid();
    bool is_OOC = aConfigurations.GetIsOOC();
    int max_rank = aConfigurations.GetMaxRank();

    // Set the floating point precision based on the template type
    FloatPoint float_point;
    if (sizeof(T) == SIZE_OF_FLOAT) {
        float_point = EXAGEOSTAT_REAL_FLOAT;
    } else if (sizeof(T) == SIZE_OF_DOUBLE) {
        float_point = EXAGEOSTAT_REAL_DOUBLE;
    } else {
        throw runtime_error("Unsupported for now!");
    }

    int MBD = lts;
    int NBD = lts;
    int MD = N;
    int ND = MBD;
    aData.GetDescriptorData()->SetDescriptor(common::HICMA_DESCRIPTOR, DESCRIPTOR_CD, is_OOC, nullptr, float_point,
                                             MBD, NBD, MBD * NBD, MD, ND, 0, 0, MD, ND, p_grid, q_grid);
    int MBUV = lts;
    int NBUV = 2 * max_rank;
    int MUV;
    int N_over_lts_times_lts = N / lts * lts;
    if (N_over_lts_times_lts < N) {
        MUV = N_over_lts_times_lts + lts;
    } else if (N_over_lts_times_lts == N) {
        MUV = N_over_lts_times_lts;
    } else {
        throw runtime_error("This case can't happens, N need to be >= lts*lts");
    }
    int expr = MUV / lts;
    int NUV = 2 * expr * max_rank;
    aData.GetDescriptorData()->SetDescriptor(common::HICMA_DESCRIPTOR, DESCRIPTOR_CUV, is_OOC, nullptr, float_point,
                                             MBUV, NBUV, MBUV * NBUV, MUV, NUV, 0, 0, MUV, NUV, p_grid, q_grid);
    auto *HICMA_descCUV = aData.GetDescriptorData()->GetDescriptor(DescriptorType::HICMA_DESCRIPTOR,
                                                                   DescriptorName::DESCRIPTOR_CUV).hicma_desc;
    int MBrk = 1;
    int NBrk = 1;
    int Mrk = HICMA_descCUV->mt;
    int Nrk = HICMA_descCUV->mt;
    aData.GetDescriptorData()->SetDescriptor(common::HICMA_DESCRIPTOR, DESCRIPTOR_CRK, is_OOC, nullptr, float_point,
                                             MBrk, NBrk, MBrk * NBrk, Mrk, Nrk, 0, 0, Mrk, Nrk, p_grid, q_grid);
}

template<typename T>
T HicmaImplementation<T>::ExaGeoStatMLETile(const hardware::ExaGeoStatHardware &aHardware, ExaGeoStatData<T> &aData,
                                            Configurations &aConfigurations, const double *theta,
                                            T *apMeasurementsMatrix) {

    this->SetContext(aHardware.GetContext(aConfigurations.GetComputation()));
    if (!aData.GetDescriptorData()->GetIsDescriptorInitiated()) {
        this->InitiateDescriptors(aConfigurations, *aData.GetDescriptorData(), apMeasurementsMatrix);
    }
    // Create a Hicma sequence, if not initialized before through the same descriptors
    RUNTIME_request_t request_array[2] = {HICMA_REQUEST_INITIALIZER, HICMA_REQUEST_INITIALIZER};
    if (!aData.GetDescriptorData()->GetSequence()) {
        HICMA_sequence_t *sequence;
        this->ExaGeoStatCreateSequence(&sequence);
        aData.GetDescriptorData()->SetSequence(sequence);
        aData.GetDescriptorData()->SetRequest(request_array);
    }
    auto pSequence = (HICMA_sequence_t *) aData.GetDescriptorData()->GetSequence();

    //Initialization
    T loglik, logdet, test_time, variance, variance1 = 1, variance2 = 1, variance3, dot_product, dot_product1, dot_product2, dot_product3, dzcpy_time, time_facto, time_solve, logdet_calculate, matrix_gen_time;
    double avg_executed_time_per_iteration = 0, avg_flops_per_iter = 0.0;

    int NRHS, i;
    T flops = 0.0;

    int N;
    int lts;
    int max_rank = aConfigurations.GetMaxRank();
    int iter_count = aData.GetMleIterations();
    auto kernel_name = aConfigurations.GetKernelName();
    int num_params = kernels::KernelsConfigurations::GetParametersNumberKernelMap()[kernel_name];
    int acc = aConfigurations.GetAccuracy();

    if (iter_count == 0) {
        this->SetModelingDescriptors(aData, aConfigurations);
    }
    auto *HICMA_descCUV = aData.GetDescriptorData()->GetDescriptor(DescriptorType::HICMA_DESCRIPTOR,
                                                                   DescriptorName::DESCRIPTOR_CUV).hicma_desc;
    auto *HICMA_descC = aData.GetDescriptorData()->GetDescriptor(DescriptorType::HICMA_DESCRIPTOR,
                                                                 DescriptorName::DESCRIPTOR_C).hicma_desc;
    auto *HICMA_descCD = aData.GetDescriptorData()->GetDescriptor(DescriptorType::HICMA_DESCRIPTOR,
                                                                  DescriptorName::DESCRIPTOR_CD).hicma_desc;
    auto *HICMA_descCrk = aData.GetDescriptorData()->GetDescriptor(DescriptorType::HICMA_DESCRIPTOR,
                                                                   DescriptorName::DESCRIPTOR_CRK).hicma_desc;
    auto *HICMA_descZ = aData.GetDescriptorData()->GetDescriptor(DescriptorType::HICMA_DESCRIPTOR,
                                                                 DescriptorName::DESCRIPTOR_Z).hicma_desc;
    auto *CHAM_descZ = aData.GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                DescriptorName::DESCRIPTOR_Z).chameleon_desc;
    auto *HICMA_descZcpy = aData.GetDescriptorData()->GetDescriptor(DescriptorType::HICMA_DESCRIPTOR,
                                                                    DescriptorName::DESCRIPTOR_Z_COPY).hicma_desc;
    auto *HICMA_desc_det = aData.GetDescriptorData()->GetDescriptor(DescriptorType::HICMA_DESCRIPTOR,
                                                                    DescriptorName::DESCRIPTOR_DETERMINANT).hicma_desc;
    auto *CHAM_desc_product = aData.GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                       DescriptorName::DESCRIPTOR_PRODUCT).chameleon_desc;

    N = HICMA_descCUV->m;
    NRHS = HICMA_descZ->n;
    lts = HICMA_descZ->mb;

    T *determinant = aData.GetDescriptorData()->GetDescriptorMatrix(HICMA_DESCRIPTOR, HICMA_desc_det);
    *determinant = 0;
    T *product = aData.GetDescriptorData()->GetDescriptorMatrix(CHAMELEON_DESCRIPTOR, CHAM_desc_product);
    *product = 0;

    string recovery_file = aConfigurations.GetRecoveryFile();
    if (recovery_file.empty() ||
        !(this->recover((char *) (recovery_file.c_str()), iter_count, (T *) theta, &loglik, num_params))) {
        if (iter_count == 0) {
            auto *z = new T[N];
            this->ExaGeoStatDesc2Lap(z, N, CHAM_descZ, EXAGEOSTAT_UPPER_LOWER);
            this->ExaGeoStatLap2Desc(z, N, HICMA_descZ, EXAGEOSTAT_UPPER_LOWER);
            delete[] z;
            //Save a copy of descZ into descZcpy for restoring each iteration
            ExaGeoStatLapackCopyTile(EXAGEOSTAT_UPPER_LOWER, HICMA_descZ, HICMA_descZcpy);
        }
    }
    //Matrix generation part.
    VERBOSE("LR:Generate New Covariance Matrix...")
    START_TIMING(matrix_gen_time);

    HICMA_problem_t hicma_problem;
    hicma_problem.theta = (double *) theta;
    hicma_problem.noise = 1e-4;
    hicma_problem.ndim = 2;

    hicma_problem.kernel_type =
            aConfigurations.GetDistanceMetric() == common::GREAT_CIRCLE_DISTANCE ? STARSH_SPATIAL_MATERN2_GCD
                                                                                 : STARSH_SPATIAL_MATERN2_SIMD;
    HICMA_zgenerate_problem(HICMA_STARSH_PROB_GEOSTAT, 'S', 0, N, lts, HICMA_descCUV->mt, HICMA_descCUV->nt,
                            &hicma_problem);
    int compress_diag = 0;
    HICMA_zgytlr_Tile(EXAGEOSTAT_LOWER, HICMA_descCUV, HICMA_descCD, HICMA_descCrk, 0, max_rank, pow(10, -1.0 * acc),
                      compress_diag, HICMA_descC);

    STOP_TIMING(matrix_gen_time);
    VERBOSE("Done.")
    //******************************
    VERBOSE("LR: re-Copy z...")
    START_TIMING(test_time);
    //re-store old Z
    this->ExaGeoStatLapackCopyTile(EXAGEOSTAT_UPPER_LOWER, HICMA_descZcpy, HICMA_descZ);
    STOP_TIMING(test_time);
    VERBOSE("Done.")

    //Calculate Cholesky Factorization (C=LL-1)
    VERBOSE("LR: Cholesky factorization of Sigma...")
    START_TIMING(time_facto);
    this->ExaGeoStatPotrfTile(EXAGEOSTAT_LOWER, HICMA_descCUV, 0, HICMA_descCD, HICMA_descCrk, max_rank, acc);

    STOP_TIMING(time_facto);
    flops = flops + flops_dpotrf(N);
    VERBOSE("Done.")

    //Calculate log(|C|) --> log(square(|L|))
    VERBOSE("LR:Calculating the log determinant ...")
    START_TIMING(logdet_calculate);
    ExaGeoStatMeasureDetTileAsync(HICMA_descCD, pSequence, &request_array[0], HICMA_desc_det);
    ExaGeoStatSequenceWait(pSequence);

    logdet = 2 * (*determinant);
    STOP_TIMING(logdet_calculate);
    VERBOSE("Done.")

    //Solving Linear System (L*X=Z)--->inv(L)*Z
    VERBOSE("LR:Solving the linear system ...")
    START_TIMING(time_solve);
    //Compute triangular solve LC*X = Z
    this->ExaGeoStatTrsmTile(EXAGEOSTAT_LEFT, EXAGEOSTAT_LOWER, EXAGEOSTAT_NO_TRANS, EXAGEOSTAT_NON_UNIT, 1,
                             HICMA_descCUV, HICMA_descCD, HICMA_descCrk, HICMA_descZ, max_rank);
    STOP_TIMING(time_solve);
    flops = flops + flops_dtrsm(ChamLeft, N, NRHS);
    VERBOSE("Done.")

    VERBOSE("LR:Calculating dot product...")
    CHAMELEON_dgemm_Tile(ChamTrans, ChamNoTrans, 1, CHAM_descZ, CHAM_descZ, 0, CHAM_desc_product);
    dot_product = *product;
    loglik = -0.5 * dot_product - 0.5 * logdet - (double) (N / 2.0) * log(2.0 * PI);
    VERBOSE("Done.")

    LOGGER(iter_count + 1 << " - Model Parameters (", true)
    if (aConfigurations.GetLogger()) {
        fprintf(aConfigurations.GetFileLogPath(), " %3d- Model Parameters (", iter_count + 1);
    }
    if ((aConfigurations.GetKernelName() == "bivariate_matern_parsimonious_profile") ||
        (aConfigurations.GetKernelName() == "bivariate_matern_parsimonious2_profile")) {
        LOGGER(setprecision(8) << variance1 << setprecision(8) << variance2)
        if (aConfigurations.GetLogger()) {
            fprintf(aConfigurations.GetFileLogPath(), "%.8f, %.8f,", variance1, variance2);
        }
        i = 2;
    } else {
        i = 0;
    }

    for (; i < num_params; i++) {
        LOGGER_PRECISION(theta[i])
        if (i < num_params - 1) {
            LOGGER_PRECISION(", ")
        }
        if (aConfigurations.GetLogger()) {
            fprintf(aConfigurations.GetFileLogPath(), "%.8f, ", theta[i]);
        }
    }
    LOGGER_PRECISION(")----> LogLi: " << loglik << "\n", 18)

    if (aConfigurations.GetLogger()) {
        fprintf(aConfigurations.GetFileLogPath(), ")----> LogLi: %.18f\n", loglik);
    }
    LOGGER(" ---- Facto Time: " << time_facto)
    LOGGER(" ---- Matrix Generation Time: " << matrix_gen_time)
    LOGGER(" ---- Total Time: " << matrix_gen_time + time_facto + logdet_calculate + time_solve)

    aData.SetMleIterations(aData.GetMleIterations() + 1);
    // for experiments
    if (Configurations::GetVerbosity() == DETAILED_MODE) {
        avg_executed_time_per_iteration += +time_facto + logdet_calculate + time_solve;
        avg_flops_per_iter += flops / 1e9 / (time_facto + time_solve);
    }

    results::Results::GetInstance()->SetMLEIterations(iter_count + 1);
    results::Results::GetInstance()->SetMaximumTheta(vector<double>(theta, theta + num_params));
    results::Results::GetInstance()->SetLogLikValue(loglik);
    return loglik;
}


template<typename T>
void HicmaImplementation<T>::ExaGeoStatLapackCopyTile(const UpperLower &aUpperLower, void *apA, void *apB) {
    int status = HICMA_dlacpy_Tile(aUpperLower, (HICMA_desc_t *) apA, (HICMA_desc_t *) apB);
    if (status != HICMA_SUCCESS) {
        throw std::runtime_error("CHAMELEON_dlacpy_Tile Failed!");
    }
}

template<typename T>
void
HicmaImplementation<T>::ExaGeoStatOptionsInit(void *apOptoins, void *apContext, void *apSequence, void *apRequest) {
    HICMA_RUNTIME_options_init((HICMA_option_t *) apOptoins, (HICMA_context_t *) apContext,
                               (HICMA_sequence_t *) apSequence, (HICMA_request_t *) apRequest);
}

template<typename T>
void HicmaImplementation<T>::ExaGeoStatOptionsFree(void *apOptions) {
    HICMA_RUNTIME_options_ws_free((HICMA_option_t *) apOptions);

}

template<typename T>
void HicmaImplementation<T>::ExaGeoStatSequenceWait(void *apSequence) {
    HICMA_Sequence_Wait((HICMA_sequence_t *) apSequence);
}

template<typename T>
void HicmaImplementation<T>::ExaGeoStatCreateSequence(void *apSequence) {
    int status = HICMA_Sequence_Create((HICMA_sequence_t **) apSequence);
    if (status != HICMA_SUCCESS) {
        throw std::runtime_error("HICMA_Sequence_Create Failed!");
    }
}

template<typename T>
void HicmaImplementation<T>::ExaGeoStatOptionsFinalize(void *apOptions, void *apContext) {
    RUNTIME_options_finalize((RUNTIME_option_t *) apOptions, (CHAM_context_t *) apContext);

}

template<typename T>
void HicmaImplementation<T>::ExaGeoStatPotrfTile(const common::UpperLower &aUpperLower, void *apA, int aDiagThick,
                                                 void *apCD, void *apCrk, const int &aMaxRank, const int &aAcc) {
    int status = HICMA_dpotrf_Tile(EXAGEOSTAT_LOWER, (HICMA_desc_t *) apA, (HICMA_desc_t *) apCD,
                                   (HICMA_desc_t *) apCrk, aDiagThick, aMaxRank, pow(10, -1.0 * aAcc));
    if (status != HICMA_SUCCESS) {
        throw std::runtime_error("HICMA_dpotrf_Tile Failed, Matrix is not positive definite");
    }

}

template<typename T>
void HicmaImplementation<T>::ExaGeoStatTrsmTile(const common::Side &aSide, const common::UpperLower &aUpperLower,
                                                const common::Trans &aTrans, const common::Diag &aDiag, const T &aAlpha,
                                                void *apA, void *apCD, void *apCrk, void *apZ, const int &aMaxRank) {
    int status = HICMA_dtrsmd_Tile(aSide, aUpperLower, aTrans, aDiag, aAlpha, (HICMA_desc_t *) apA,
                                   (HICMA_desc_t *) apCD, (HICMA_desc_t *) apCrk, (HICMA_desc_t *) apZ, aMaxRank);
    if (status != HICMA_SUCCESS) {
        throw std::runtime_error("HICMA_dtrsmd_Tile Failed!");
    }
}

template<typename T>
int HicmaImplementation<T>::ExaGeoStatMeasureDetTileAsync(void *apDescA, void *apSequence, void *apRequest,
                                                          void *apDescDet) {

    // Check for initialize the Hicma context.
    if (!this->mpContext) {
        throw std::runtime_error(
                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStatHardware(computation, cores_number, gpu_numbers);'.");
    }
    HICMA_option_t options;
    this->ExaGeoStatOptionsInit(&options, this->mpContext, apSequence, apRequest);

    int m;
    int tempmm;
    auto Z = (HICMA_desc_t *) apDescA;
    auto det = (HICMA_desc_t *) apDescDet;
    struct starpu_codelet *cl = &this->cl_dmdet;

    for (m = 0; m < Z->mt; m++) {
        tempmm = m == Z->mt - 1 ? Z->m - m * Z->mb : Z->mb;
        starpu_insert_task(cl,
                           STARPU_VALUE, &tempmm, sizeof(int),
                           STARPU_VALUE, &tempmm, sizeof(int),
                           STARPU_R, ExaGeoStatDataGetAddr(Z, m, 0),
                           STARPU_RW, ExaGeoStatDataGetAddr(det, 0, 0),
                           0);
    }
    this->ExaGeoStatOptionsFree(&options);
    this->ExaGeoStatOptionsFinalize(&options, (HICMA_context_t * )
    this->mpContext);
    return HICMA_SUCCESS;
}


template<typename T>
void HicmaImplementation<T>::ExaGeoStatLap2Desc(T *apA, const int &aLDA, void *apDescA, const UpperLower &aUpperLower) {
    int status = HICMA_Lapack_to_Tile(apA, aLDA, (HICMA_desc_t *) apDescA);
    if (status != HICMA_SUCCESS) {
        throw std::runtime_error("HICMA_Lapack_to_Tile Failed!");
    }
}

template<typename T>
void *HicmaImplementation<T>::ExaGeoStatDataGetAddr(void *apA, int aAm, int aAn) {
    return HICMA_RUNTIME_data_getaddr((HICMA_desc_t *) apA, aAm, aAn);

}