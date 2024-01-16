
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file HicmaImplementation.cpp
 * @brief Sets up the HiCMA descriptors needed for the tile low rank computations in ExaGeoStat.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-03-26
**/

#include <linear-algebra-solvers/concrete/hicma/tlr/HicmaImplementation.hpp>

using namespace std;

using namespace exageostat::linearAlgebra::tileLowRank;
using namespace exageostat::common;
using namespace exageostat::dataunits;
using namespace exageostat::kernels;
using namespace exageostat::helpers;
using namespace exageostat::hardware;
using namespace exageostat::configurations;

int store_only_diagonal_tiles = 1;
int use_scratch = 1;
int global_check = 0;  //used to create dense matrix for accuracy check

template<typename T>
void HicmaImplementation<T>::SetModelingDescriptors(std::unique_ptr<dataunits::ExaGeoStatData<T>> &aData,
                                                    Configurations &aConfigurations, const int &aP) {

    int full_problem_size = aConfigurations.GetProblemSize() * aP;
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
    int MD = full_problem_size;
    int ND = MBD;
    aData->GetDescriptorData()->SetDescriptor(common::HICMA_DESCRIPTOR, DESCRIPTOR_CD, is_OOC, nullptr, float_point,
                                              MBD,
                                              NBD, MBD * NBD, MD, ND, 0, 0, MD, ND, p_grid, q_grid);
    int MBUV = lts;
    int NBUV = 2 * max_rank;
    int MUV;
    int N_over_lts_times_lts = full_problem_size / lts * lts;
    if (N_over_lts_times_lts < full_problem_size) {
        MUV = N_over_lts_times_lts + lts;
    } else if (N_over_lts_times_lts == full_problem_size) {
        MUV = N_over_lts_times_lts;
    } else {
        throw runtime_error("This case can't happens, N need to be >= lts*lts");
    }
    int expr = MUV / lts;
    int NUV = 2 * expr * max_rank;
    aData->GetDescriptorData()->SetDescriptor(common::HICMA_DESCRIPTOR, DESCRIPTOR_CUV, is_OOC, nullptr, float_point,
                                              MBUV, NBUV, MBUV * NBUV, MUV, NUV, 0, 0, MUV, NUV, p_grid, q_grid);
    auto *HICMA_descCUV = aData->GetDescriptorData()->GetDescriptor(DescriptorType::HICMA_DESCRIPTOR,
                                                                    DescriptorName::DESCRIPTOR_CUV).hicma_desc;
    int MBrk = 1;
    int NBrk = 1;
    int Mrk = HICMA_descCUV->mt;
    int Nrk = HICMA_descCUV->mt;
    aData->GetDescriptorData()->SetDescriptor(common::HICMA_DESCRIPTOR, DESCRIPTOR_CRK, is_OOC, nullptr, float_point,
                                              MBrk, NBrk, MBrk * NBrk, Mrk, Nrk, 0, 0, Mrk, Nrk, p_grid, q_grid);

    aData->GetDescriptorData()->SetDescriptor(common::HICMA_DESCRIPTOR, DESCRIPTOR_Z, is_OOC, nullptr, float_point,
                                              lts, lts, lts * lts, full_problem_size, 1, 0, 0, full_problem_size, 1, p_grid, q_grid);

    aData->GetDescriptorData()->SetDescriptor(common::HICMA_DESCRIPTOR, DESCRIPTOR_Z_COPY, is_OOC, nullptr, float_point,
                                              lts, lts, lts * lts, full_problem_size, 1, 0, 0, full_problem_size, 1, p_grid, q_grid);
    if (aConfigurations.GetIsNonGaussian()) {
        aData->GetDescriptorData()->SetDescriptor(common::HICMA_DESCRIPTOR, DESCRIPTOR_PRODUCT, is_OOC, nullptr,
                                                  float_point,
                                                  lts, lts, lts * lts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);
        aData->GetDescriptorData()->SetDescriptor(common::HICMA_DESCRIPTOR, DESCRIPTOR_SUM, is_OOC, nullptr,
                                                  float_point,
                                                  lts, lts, lts * lts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);
    }
}

template<typename T>
T HicmaImplementation<T>::ExaGeoStatMLETile(const hardware::ExaGeoStatHardware &aHardware,
                                            std::unique_ptr<dataunits::ExaGeoStatData<T>> &aData,
                                            Configurations &aConfigurations, const double *theta,
                                            T *apMeasurementsMatrix, const Kernel<T> &aKernel) {

    this->SetContext(aHardware.GetContext(aConfigurations.GetComputation()));
    if (!aData->GetDescriptorData()->GetIsDescriptorInitiated()) {
        this->InitiateDescriptors(aConfigurations, *aData->GetDescriptorData(),aKernel.GetVariablesNumber(), apMeasurementsMatrix);
    }
    // Create a Hicma sequence, if not initialized before through the same descriptors
    RUNTIME_request_t request_array[2] = {HICMA_REQUEST_INITIALIZER, HICMA_REQUEST_INITIALIZER};
    if (!aData->GetDescriptorData()->GetSequence()) {
        HICMA_sequence_t *sequence;
        this->ExaGeoStatCreateSequence(&sequence);
        aData->GetDescriptorData()->SetSequence(sequence);
        aData->GetDescriptorData()->SetRequest(request_array);
    }
    auto pSequence = (HICMA_sequence_t *) aData->GetDescriptorData()->GetSequence();

    //Initialization
    T loglik, logdet, test_time, variance, variance1 = 1, variance2 = 1, variance3, dot_product, dot_product1, dot_product2, dot_product3, dzcpy_time, time_facto, time_solve, logdet_calculate, matrix_gen_time;
    double accumulated_executed_time, accumulated_flops;

    int NRHS, i;
    T flops = 0.0;

    int N;
    int lts;
    int max_rank = aConfigurations.GetMaxRank();
    int iter_count = aData->GetMleIterations();
    auto kernel_name = aConfigurations.GetKernelName();
    int num_params = aKernel.GetParametersNumbers();
    int acc = aConfigurations.GetAccuracy();

    if (iter_count == 0) {
        this->SetModelingDescriptors(aData, aConfigurations, aKernel.GetVariablesNumber());
    }
    auto *HICMA_descCUV = aData->GetDescriptorData()->GetDescriptor(DescriptorType::HICMA_DESCRIPTOR,
                                                                    DescriptorName::DESCRIPTOR_CUV).hicma_desc;
    auto *HICMA_descC = aData->GetDescriptorData()->GetDescriptor(DescriptorType::HICMA_DESCRIPTOR,
                                                                  DescriptorName::DESCRIPTOR_C).hicma_desc;
    auto *HICMA_descCD = aData->GetDescriptorData()->GetDescriptor(DescriptorType::HICMA_DESCRIPTOR,
                                                                   DescriptorName::DESCRIPTOR_CD).hicma_desc;
    auto *HICMA_descCrk = aData->GetDescriptorData()->GetDescriptor(DescriptorType::HICMA_DESCRIPTOR,
                                                                    DescriptorName::DESCRIPTOR_CRK).hicma_desc;
    auto *HICMA_descZ = aData->GetDescriptorData()->GetDescriptor(DescriptorType::HICMA_DESCRIPTOR,
                                                                  DescriptorName::DESCRIPTOR_Z).hicma_desc;
    auto *CHAM_descZ = aData->GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                 DescriptorName::DESCRIPTOR_Z).chameleon_desc;
    auto *CHAM_descZcpy = aData->GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                    DescriptorName::DESCRIPTOR_Z_COPY).chameleon_desc;
    auto *HICMA_descZcpy = aData->GetDescriptorData()->GetDescriptor(DescriptorType::HICMA_DESCRIPTOR,
                                                                     DescriptorName::DESCRIPTOR_Z_COPY).hicma_desc;
    auto *HICMA_desc_det = aData->GetDescriptorData()->GetDescriptor(DescriptorType::HICMA_DESCRIPTOR,
                                                                     DescriptorName::DESCRIPTOR_DETERMINANT).hicma_desc;
    auto *CHAM_desc_product = aData->GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                        DescriptorName::DESCRIPTOR_PRODUCT).chameleon_desc;
    auto *HICMA_desc_product = aData->GetDescriptorData()->GetDescriptor(DescriptorType::HICMA_DESCRIPTOR,
                                                                         DescriptorName::DESCRIPTOR_PRODUCT).chameleon_desc;
    auto *HICMA_desc_sum = aData->GetDescriptorData()->GetDescriptor(DescriptorType::HICMA_DESCRIPTOR,
                                                                     DescriptorName::DESCRIPTOR_SUM).chameleon_desc;
    N = HICMA_descCUV->m;
    NRHS = HICMA_descZ->n;
    lts = HICMA_descZ->mb;

    T *determinant = aData->GetDescriptorData()->GetDescriptorMatrix(HICMA_DESCRIPTOR, HICMA_desc_det);
    *determinant = 0;
    T *product = aData->GetDescriptorData()->GetDescriptorMatrix(CHAMELEON_DESCRIPTOR, CHAM_desc_product);
    *product = 0;
    T *sum;

    string recovery_file = aConfigurations.GetRecoveryFile();
    if (recovery_file.empty() ||
        !(this->recover((char *) (recovery_file.c_str()), iter_count, (T *) theta, &loglik, num_params))) {
        if (iter_count == 0) {
            // Copy dense matrix Z into tile low rank Z.
            this->CopyDescriptors(CHAM_descZ, HICMA_descZ, N, CHAMELEON_TO_HICMA);
            // Save a copy of descZ into descZcpy for restoring each iteration
            this->ExaGeoStatLapackCopyTile(EXAGEOSTAT_UPPER_LOWER, HICMA_descZ, HICMA_descZcpy);
            // Save another copy into descZcpy for chameleon, This is in case of other operations after Modeling. ex: Prediction.
            CHAMELEON_dlacpy_Tile(ChamUpperLower, CHAM_descZ, CHAM_descZcpy);
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
    int hicma_data_type;
    if (aConfigurations.GetIsNonGaussian()) {
        hicma_data_type = HICMA_STARSH_PROB_GEOSTAT_NON_GAUSSIAN;
    } else {
        hicma_data_type = HICMA_STARSH_PROB_GEOSTAT;
    }

    HICMA_zgenerate_problem(hicma_data_type, 'S', 0, N, lts, HICMA_descCUV->mt, HICMA_descCUV->nt,
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
    this->ExaGeoStatPotrfTile(EXAGEOSTAT_LOWER, HICMA_descCUV, 0, HICMA_descCD, HICMA_descCrk, max_rank,
                              pow(10, -1.0 * acc));

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

    if (aConfigurations.GetIsNonGaussian()) {
        VERBOSE("Transform Z vector to Gaussian field ...")
        this->ExaGeoStatNonGaussianTransformTileAsync(HICMA_descZ, HICMA_desc_product, (T *) theta, pSequence,
                                                      &request_array[0]);
        ExaGeoStatSequenceWait(pSequence);

        VERBOSE(" Done.")
        sum = aData->GetDescriptorData()->GetDescriptorMatrix(HICMA_DESCRIPTOR, HICMA_desc_sum);
        *sum = 0;

        VERBOSE("LR:Calculating dot product...")
        CHAMELEON_dgemm_Tile(ChamTrans, ChamNoTrans, 1, CHAM_descZ, CHAM_descZ, 0, CHAM_desc_product);
        VERBOSE(" Done.")

        double global_sum;
#if defined(CHAMELEON_USE_MPI)
        MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
        global_sum = *product;
#endif
        this->ExaGeoStatNonGaussianLogLikeTileAsync(HICMA_descZ, HICMA_desc_sum, (T *) theta, pSequence,
                                                    &request_array[0]);
        ExaGeoStatSequenceWait(pSequence);
        VERBOSE(" Done.")
    }

    //Solving Linear System (L*X=Z)--->inv(L)*Z
    VERBOSE("LR:Solving the linear system ...")
    START_TIMING(time_solve);
    //Compute triangular solve LC*X = Z
    this->ExaGeoStatTrsmTile(EXAGEOSTAT_LEFT, EXAGEOSTAT_LOWER, EXAGEOSTAT_NO_TRANS, EXAGEOSTAT_NON_UNIT, 1,
                             HICMA_descCUV, HICMA_descCD, HICMA_descCrk, HICMA_descZ, max_rank);
    STOP_TIMING(time_solve);
    flops = flops + flops_dtrsm(ChamLeft, N, NRHS);
    VERBOSE("Done.")

    VERBOSE("Copy to chameleon")
    this->CopyDescriptors(HICMA_descZ, CHAM_descZ, N, HICMA_TO_CHAMELEON);

    VERBOSE("LR:Calculating dot product...")
    CHAMELEON_dgemm_Tile(ChamTrans, ChamNoTrans, 1, CHAM_descZ, CHAM_descZ, 0, CHAM_desc_product);
    dot_product = *product;
    loglik = -0.5 * dot_product - 0.5 * logdet;
    if (aConfigurations.GetIsNonGaussian()) {
        loglik = loglik - *sum - N * log(theta[3]) - (double) (N / 2.0) * log(2.0 * PI);
    } else {
        loglik = loglik - (double) (N / 2.0) * log(2.0 * PI);
    }
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
    LOGGER(" ---- Log Determent Time: " << logdet_calculate)
    LOGGER(" ---- dtrsm Time: " << time_solve)
    LOGGER(" ---- Matrix Generation Time: " << matrix_gen_time)
    LOGGER(" ---- Total Time: " << time_facto + logdet_calculate + time_solve)
    LOGGER(" ---- Gflop/s: " << flops / 1e9 / (time_facto + time_solve))

    aData->SetMleIterations(aData->GetMleIterations() + 1);

    // for experiments and benchmarking
    accumulated_executed_time =
            results::Results::GetInstance()->GetTotalModelingExecutionTime() + time_facto + logdet_calculate +
            time_solve;
    results::Results::GetInstance()->SetTotalModelingExecutionTime(accumulated_executed_time);
    accumulated_flops =
            results::Results::GetInstance()->GetTotalModelingFlops() + (flops / 1e9 / (time_facto + time_solve));
    results::Results::GetInstance()->SetTotalModelingFlops(accumulated_flops);

    results::Results::GetInstance()->SetMLEIterations(iter_count + 1);
    results::Results::GetInstance()->SetMaximumTheta(vector<double>(theta, theta + num_params));
    results::Results::GetInstance()->SetLogLikValue(loglik);

    aConfigurations.SetEstimatedTheta(aConfigurations.GetStartingTheta());
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
HicmaImplementation<T>::ExaGeoStatOptionsInit(void *apOptions, void *apContext, void *apSequence, void *apRequest) {
    HICMA_RUNTIME_options_init((HICMA_option_t *) apOptions, (HICMA_context_t *) apContext,
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
void HicmaImplementation<T>::ExaGeoStatPotrfTile(const common::UpperLower &aUpperLower, void *apA, int aBand,
                                                 void *apCD, void *apCrk, const int &aMaxRank, const int &aAcc) {
    int status = HICMA_dpotrf_Tile(EXAGEOSTAT_LOWER, (HICMA_desc_t *) apA, (HICMA_desc_t *) apCD,
                                   (HICMA_desc_t *) apCrk, aBand, aMaxRank, pow(10, -1.0 * aAcc));
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
    int temp;
    auto Z = (HICMA_desc_t *) apDescA;
    auto det = (HICMA_desc_t *) apDescDet;
    struct starpu_codelet *cl = &this->cl_dmdet;

    for (m = 0; m < Z->mt; m++) {
        temp = m == Z->mt - 1 ? Z->m - m * Z->mb : Z->mb;
        starpu_insert_task(cl,
                           STARPU_VALUE, &temp, sizeof(int),
                           STARPU_VALUE, &temp, sizeof(int),
                           STARPU_R, ExaGeoStatDataGetAddr(Z, m, 0),
                           STARPU_RW, ExaGeoStatDataGetAddr(det, 0, 0),
                           0);
    }
    this->ExaGeoStatOptionsFree(&options);
    this->ExaGeoStatOptionsFinalize(&options, (HICMA_context_t *)
            this->mpContext);
    return HICMA_SUCCESS;
}

template<typename T>
void *HicmaImplementation<T>::ExaGeoStatDataGetAddr(void *apA, int aAm, int aAn) {
    return HICMA_RUNTIME_data_getaddr((HICMA_desc_t *) apA, aAm, aAn);

}

template<typename T>
int
HicmaImplementation<T>::ExaGeoStatNonGaussianLogLikeTileAsync(void *apDescZ, void *apDescSum, const T *apTheta,
                                                              void *apSequence,
                                                              void *apRequest) {
    // Check for initialize the Hicma context.
    if (!this->mpContext) {
        throw std::runtime_error(
                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStatHardware(computation, cores_number, gpu_numbers);'.");
    }
    HICMA_option_t options;
    this->ExaGeoStatOptionsInit(&options, this->mpContext, apSequence, apRequest);

    int m, m0;
    int tempmm;
    auto Z = (HICMA_desc_t *) apDescZ;

    struct starpu_codelet *cl = &this->cl_non_gaussian_loglike_lr;


    for (m = 0; m < Z->mt; m++) {
        tempmm = m == Z->mt - 1 ? Z->m - m * Z->mb : Z->mb;

        m0 = m * Z->mb;

        starpu_insert_task(cl,
                           STARPU_VALUE, &tempmm, sizeof(int),
                           STARPU_VALUE, &m0, sizeof(int),
                           STARPU_R, ExaGeoStatDataGetAddr(apDescZ, m, 0),
                           STARPU_RW, ExaGeoStatDataGetAddr(apDescSum, 0, 0),
                           STARPU_VALUE, &apTheta[0], sizeof(double),
                           STARPU_VALUE, &apTheta[1], sizeof(double),
                           STARPU_VALUE, &apTheta[2], sizeof(double),
                           STARPU_VALUE, &apTheta[3], sizeof(double),
                           STARPU_VALUE, &apTheta[4], sizeof(double),
                           STARPU_VALUE, &apTheta[5], sizeof(double),
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
                STARPU_NAME, "non_gaussian_loglike_lr",
#endif
                           0);
    }

    this->ExaGeoStatOptionsFree(&options);
    this->ExaGeoStatOptionsFinalize(&options, (HICMA_context_t *)
            this->mpContext);
    return HICMA_SUCCESS;
}


template<typename T>
int HicmaImplementation<T>::ExaGeoStatNonGaussianTransformTileAsync(void *apDescZ, void *apDescFlag, const T *apTheta,
                                                                    void *apSequence, void *apRequest) {
    // Check for initialize the Hicma context.
    if (!this->mpContext) {
        throw std::runtime_error(
                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStatHardware(computation, cores_number, gpu_numbers);'.");
    }
    HICMA_option_t options;
    this->ExaGeoStatOptionsInit(&options, this->mpContext, apSequence, apRequest);

    int m, m0;
    int tempmm;
    auto Z = (HICMA_desc_t *) apDescZ;
    struct starpu_codelet *cl = &this->cl_non_gaussian_transform_lr;


    for (m = 0; m < Z->mt; m++) {
        tempmm = m == Z->mt - 1 ? Z->m - m * Z->mb : Z->mb;
        m0 = m * Z->mb;
        starpu_insert_task(cl,
                           STARPU_VALUE, &tempmm, sizeof(int),
                           STARPU_VALUE, &m0, sizeof(int),
                           STARPU_RW, ExaGeoStatDataGetAddr(apDescZ, m, 0),
                           STARPU_VALUE, &apTheta[0], sizeof(double),
                           STARPU_VALUE, &apTheta[1], sizeof(double),
                           STARPU_VALUE, &apTheta[2], sizeof(double),
                           STARPU_VALUE, &apTheta[3], sizeof(double),
                           STARPU_VALUE, &apTheta[4], sizeof(double),
                           STARPU_VALUE, &apTheta[5], sizeof(double),
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
                STARPU_NAME, "non_gaussian_transform_lr",
#endif
                           0);
    }
    this->ExaGeoStatOptionsFree(&options);
    this->ExaGeoStatOptionsFinalize(&options, (HICMA_context_t *)
            this->mpContext);
    return HICMA_SUCCESS;
}