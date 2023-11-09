
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ChameleonImplementation.cpp
 * @brief This file contains the declaration of ChameleonImplementation class.
 * @details ChameleonImplementation is a concrete implementation of LinearAlgebraMethods class for dense or diagonal-super tile matrices..
 * @version 1.0.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-03-20
**/

#ifdef USE_MPI
#include <mpi.h>
#endif

#include <cblas.h>
#include <linear-algebra-solvers/concrete/chameleon/ChameleonImplementation.hpp>

using namespace std;

using namespace exageostat::common;
using namespace exageostat::dataunits;
using namespace exageostat::configurations;
using namespace exageostat::linearAlgebra;

template<typename T>
int ChameleonImplementation<T>::ExaGeoStatDoubleDotProduct(void *apDescA, void *apDescProduct, void *apSequence,
                                                            void *apRequest) {
    RUNTIME_option_t options;
    this->ExaGeoStatOptionsInit(&options, this->mpContext, apSequence, apRequest);


    int m, m0;
    int tempmm;
    auto A = (CHAM_desc_t *) apDescA;

    struct starpu_codelet *cl=&this->cl_ddotp;

    for (m = 0; m < A->mt; m++) {
        tempmm = m == A->mt-1 ? A->m - m * A->mb : A->mb;

        m0 = m * A->mb;


        starpu_insert_task(cl,
                           STARPU_VALUE, &tempmm, sizeof(int),
                           STARPU_VALUE, &m0,   sizeof(int),
                           STARPU_RW, ExaGeoStatDataGetAddr(apDescProduct, 0, 0),
                           STARPU_R, RUNTIME_data_getaddr(A, m, 0),
                           0);
    }
    this->ExaGeoStatOptionsFree(&options);
    this->ExaGeoStatOptionsFinalize(&options, (CHAM_context_t *) this->mpContext);
    return CHAMELEON_SUCCESS;
}


template<typename T>
T ChameleonImplementation<T>::ExaGeoStatMLETile(const hardware::ExaGeoStatHardware &aHardware, ExaGeoStatData<T> &aData,
                                                Configurations &aConfigurations, const double *theta,
                                                T *apMeasurementsMatrix, const kernels::Kernel<T> &aKernel) {

    this->SetContext(aHardware.GetContext(aConfigurations.GetComputation()));
    if (!aData.GetDescriptorData()->GetIsDescriptorInitiated()) {
        this->InitiateDescriptors(aConfigurations, *aData.GetDescriptorData(), apMeasurementsMatrix);
    }
    // Create a Chameleon sequence, if not initialized before through the same descriptors
    RUNTIME_request_t request_array[2] = {RUNTIME_REQUEST_INITIALIZER, RUNTIME_REQUEST_INITIALIZER};
    if (!aData.GetDescriptorData()->GetSequence()) {
        RUNTIME_sequence_t *sequence;
        this->ExaGeoStatCreateSequence(&sequence);
        aData.GetDescriptorData()->SetSequence(sequence);
        aData.GetDescriptorData()->SetRequest(request_array);
    }

    auto pSequence = (RUNTIME_sequence_t *) aData.GetDescriptorData()->GetSequence();
    //Initialization
    T loglik = 0.0, logdet, variance, variance1 = 1, variance2 = 1, variance3, dot_product = 0, dot_product1 = 0, dot_product2 = 0, dot_product3, n, dzcpy_time, time_facto, time_solve, logdet_calculate, matrix_gen_time;
    double accumulated_executed_time, accumulated_flops;

    int nhrs, i;
    T flops = 0.0;
    T *univariate_theta, *univariate2_theta, *univariate3_theta, nu12, rho, sigma_square12;

    auto kernel_name = aConfigurations.GetKernelName();
    int num_params = aKernel.GetParametersNumbers();
    auto median_locations = Locations<T>(1, aData.GetLocations()->GetDimension());
    aData.CalculateMedianLocations(kernel_name, median_locations);

    auto *CHAM_desc_C = aData.GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                 DescriptorName::DESCRIPTOR_C).chameleon_desc;
    auto *CHAM_desc_sub_C11 = aData.GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                       DescriptorName::DESCRIPTOR_C11).chameleon_desc;
    auto *CHAM_desc_sub_C12 = aData.GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                       DescriptorName::DESCRIPTOR_C12).chameleon_desc;
    auto *CHAM_desc_sub_C22 = aData.GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                       DescriptorName::DESCRIPTOR_C22).chameleon_desc;
    auto *CHAM_desc_Z = aData.GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                 DescriptorName::DESCRIPTOR_Z).chameleon_desc;
    auto *CHAM_desc_Z1 = aData.GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                  DescriptorName::DESCRIPTOR_Z_1).chameleon_desc;
    auto *CHAM_desc_Z2 = aData.GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                  DescriptorName::DESCRIPTOR_Z_2).chameleon_desc;
    auto *CHAM_desc_Z3 = aData.GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                  DescriptorName::DESCRIPTOR_Z_3).chameleon_desc;
    auto *CHAM_desc_Zcpy = aData.GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                    DescriptorName::DESCRIPTOR_Z_COPY).chameleon_desc;
    auto *CHAM_desc_det = aData.GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                   DescriptorName::DESCRIPTOR_DETERMINANT).chameleon_desc;
    auto *CHAM_desc_product = aData.GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                       DescriptorName::DESCRIPTOR_PRODUCT).chameleon_desc;
    auto *CHAM_desc_product1 = aData.GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                        DescriptorName::DESCRIPTOR_PRODUCT_1).chameleon_desc;
    auto *CHAM_desc_product2 = aData.GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                        DescriptorName::DESCRIPTOR_PRODUCT_2).chameleon_desc;
    auto *CHAM_desc_product3 = aData.GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                        DescriptorName::DESCRIPTOR_PRODUCT_3).chameleon_desc;

    T *determinant = aData.GetDescriptorData()->GetDescriptorMatrix(CHAMELEON_DESCRIPTOR, CHAM_desc_det);
    *determinant = 0;
    T *product = aData.GetDescriptorData()->GetDescriptorMatrix(CHAMELEON_DESCRIPTOR, CHAM_desc_product);
    *product = 0;
    n = CHAM_desc_C->m;
    nhrs = CHAM_desc_Z->n;

    string recovery_file = aConfigurations.GetRecoveryFile();
    int iter_count = aData.GetMleIterations();

    if (recovery_file.empty() ||
        !(this->recover((char *) (recovery_file.c_str()), iter_count, (T *) theta, &loglik, num_params))) {
        START_TIMING(dzcpy_time);
        if (iter_count == 0) {
            // Save a copy of descZ into descZcpy for restoring each iteration (Only for the first iteration)
            this->ExaGeoStatLapackCopyTile(EXAGEOSTAT_UPPER_LOWER, CHAM_desc_Z, CHAM_desc_Zcpy);
        } else {
            VERBOSE("Re-store the original Z vector...")
            this->ExaGeoStatLapackCopyTile(EXAGEOSTAT_UPPER_LOWER, CHAM_desc_Zcpy, CHAM_desc_Z);
            VERBOSE("Done.")
        }
        STOP_TIMING(dzcpy_time);
    }

    //Generate new co-variance matrix C based on new theta
    VERBOSE("Generate New Covariance Matrix...")
    START_TIMING(matrix_gen_time);

    if (kernel_name == "bivariate_matern_parsimonious2" ||
        kernel_name == "bivariate_matern_parsimonious2_profile") {

        int upper_lower = EXAGEOSTAT_UPPER_LOWER;
        univariate_theta = new T[3];
        univariate2_theta = new T[3];
        univariate3_theta = new T[3];
        univariate_theta[0] = theta[0];
        univariate_theta[1] = theta[2];
        univariate_theta[2] = theta[3];

        int distance_metric = aConfigurations.GetDistanceMetric();

        this->CovarianceMatrixCodelet(*aData.GetDescriptorData(), CHAM_desc_sub_C11, upper_lower, aData.GetLocations(),
                                      aData.GetLocations(), &median_locations, univariate_theta, distance_metric,
                                      &aKernel);

        nu12 = 0.5 * (theta[3] + theta[4]);
        rho = theta[5] * sqrt((tgamma(theta[3] + 1) * tgamma(theta[4] + 1)) / (tgamma(theta[3]) * tgamma(theta[4]))) *
              tgamma(nu12) / tgamma(nu12 + 1);
        sigma_square12 = rho * sqrt(theta[0] * theta[1]);
        univariate2_theta[0] = sigma_square12;
        univariate2_theta[1] = theta[2];
        univariate2_theta[2] = nu12;

        this->CovarianceMatrixCodelet(*aData.GetDescriptorData(), CHAM_desc_sub_C12, upper_lower, &median_locations,
                                      aData.GetLocations(), &median_locations, univariate2_theta, 0, &aKernel);
        STOP_TIMING(matrix_gen_time);
        VERBOSE("Done.")

        univariate3_theta[0] = theta[1];
        univariate3_theta[1] = theta[2];
        univariate3_theta[2] = theta[4];

        this->CovarianceMatrixCodelet(*aData.GetDescriptorData(), CHAM_desc_sub_C22, upper_lower, &median_locations,
                                      aData.GetLocations(), &median_locations, univariate2_theta, 0, &aKernel);
    } else {
        int upper_lower = EXAGEOSTAT_LOWER;
        this->CovarianceMatrixCodelet(*aData.GetDescriptorData(), CHAM_desc_C, upper_lower, aData.GetLocations(),
                                      aData.GetLocations(), &median_locations, (T *) theta, 0, &aKernel);
    }
    this->ExaGeoStatSequenceWait(pSequence);
    STOP_TIMING(matrix_gen_time);
    VERBOSE("Done.")

    VERBOSE("Cholesky factorization of Sigma...")
    START_TIMING(time_facto);
    this->ExaGeoStatPotrfTile(EXAGEOSTAT_LOWER, CHAM_desc_C, aConfigurations.GetBand(), nullptr, nullptr, 0, 0);
    STOP_TIMING(time_facto);
    flops += flops_dpotrf(n);
    VERBOSE("Done.")

    //Calculate log(|C|) --> log(square(|L|))
    VERBOSE("Calculating the log determinant ...")
    START_TIMING(logdet_calculate);
    this->ExaGeoStatMeasureDetTileAsync(CHAM_desc_C, pSequence, &request_array, CHAM_desc_det);
    this->ExaGeoStatSequenceWait(pSequence);

    logdet = 2 * (*determinant);
    STOP_TIMING(logdet_calculate);
    VERBOSE("Done.")

    // Solving Linear System (L*X=Z)--->inv(L)*Z
    VERBOSE("Solving the linear system ...")
    START_TIMING(time_solve);
    this->ExaGeoStatTrsmTile(EXAGEOSTAT_LEFT, EXAGEOSTAT_LOWER, EXAGEOSTAT_NO_TRANS, EXAGEOSTAT_NON_UNIT, 1,
                             CHAM_desc_C, nullptr, nullptr, CHAM_desc_Z, 0);
    STOP_TIMING(time_solve);
    flops += flops_dtrsm(ChamLeft, n, nhrs);
    VERBOSE("Done.")

    //Calculate MLE likelihood
    VERBOSE("Calculating the MLE likelihood function ...")
    ExaGeoStatDoubleDotProduct(CHAM_desc_Z, CHAM_desc_product, pSequence, request_array);
    ExaGeoStatSequenceWait(pSequence);
    if (kernel_name == "BivariateMaternParsimonious2Profile") {
        loglik =
                -(n / 2) + (n / 2) * log(n) - (n / 2) * log(dot_product) - 0.5 * logdet - (T) (n / 2.0) * log(2.0 * PI);
        CHAMELEON_dgemm_Tile(ChamTrans, ChamNoTrans, 1, CHAM_desc_Z1, CHAM_desc_Z1, 0, CHAM_desc_product1);
        CHAMELEON_dgemm_Tile(ChamTrans, ChamNoTrans, 1, CHAM_desc_Z2, CHAM_desc_Z2, 0, CHAM_desc_product2);
        variance1 = (1.0 / (n / 2)) * dot_product1;
        variance2 = (1.0 / (n / 2)) * dot_product2;

    } else if (kernel_name == "BivariateMaternParsimoniousProfile") {
        loglik = -(n / 2) + (n / 2) * log(n) - (n / 2) * log(dot_product) - 0.5 * logdet -
                 (double) (n / 2.0) * log(2.0 * PI);
        this->ExaGeoStaStrideVectorTileAsync(CHAM_desc_Z, CHAM_desc_Z1, CHAM_desc_Z2, pSequence, &request_array[0]);
        CHAMELEON_dgemm_Tile(ChamTrans, ChamNoTrans, 1, CHAM_desc_Z1, CHAM_desc_Z1, 0, CHAM_desc_product1);
        CHAMELEON_dgemm_Tile(ChamTrans, ChamNoTrans, 1, CHAM_desc_Z2, CHAM_desc_Z2, 0, CHAM_desc_product2);
        variance1 = (1.0 / (n / 2)) * dot_product1;
        variance2 = (1.0 / (n / 2)) * dot_product2;
    } else if (kernel_name == "TrivariateMaternParsimoniousProfile") {

        loglik = -(n / 3.0) + (n / 3.0) * log(n / 3.0) - (n / 3.0) * log(dot_product) - 0.5 * logdet -
                 (double) (n / 3.0) * log(2.0 * PI);
        //to be optimized
        this->ExaGeoStaStrideVectorTileAsync(CHAM_desc_Z, CHAM_desc_Z1, CHAM_desc_Z2, CHAM_desc_Z3, pSequence,
                                             &request_array[0]);
        CHAMELEON_dgemm_Tile(ChamTrans, ChamNoTrans, 1, CHAM_desc_Z1, CHAM_desc_Z1, 0, CHAM_desc_product1);
        CHAMELEON_dgemm_Tile(ChamTrans, ChamNoTrans, 1, CHAM_desc_Z2, CHAM_desc_Z2, 0, CHAM_desc_product2);
        CHAMELEON_dgemm_Tile(ChamTrans, ChamNoTrans, 1, CHAM_desc_Z3, CHAM_desc_Z3, 0, CHAM_desc_product3);
        variance1 = (1.0 / (n / 3.0)) * dot_product1;
        variance2 = (1.0 / (n / 3.0)) * dot_product2;
    } else {
        dot_product = *product;
        loglik = -0.5 * dot_product - 0.5 * logdet - (double) (n / 2.0) * log(2.0 * PI);
    }
    VERBOSE("Done.")

    //Distribute the values in the case of MPI
#if defined(CHAMELEON_USE_MPI)
    MPI_Bcast(&loglik, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
#endif

    LOGGER(iter_count + 1 << " - Model Parameters (", true)

    if (aConfigurations.GetLogger()) {
        fprintf(aConfigurations.GetFileLogPath(), " %3d- Model Parameters (", iter_count + 1);
    }

    if ((aConfigurations.GetKernelName() == "BivariateMaternParsimoniousProfile") ||
        (aConfigurations.GetKernelName() == "BivariateMaternParsimonious2Profile")) {
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
    LOGGER(" ---- Gflop/s: " << flops / 1e9 / (time_facto  + time_solve))

    aData.SetMleIterations(aData.GetMleIterations() + 1);

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
void
ChameleonImplementation<T>::ExaGeoStatLap2Desc(T *apA, const int &aLDA, void *apDescA, const UpperLower &aUpperLower) {
    int status = CHAMELEON_Lap2Desc((cham_uplo_t) aUpperLower, apA, aLDA, (CHAM_desc_t *) apDescA);
    if (status != CHAMELEON_SUCCESS) {
        throw std::runtime_error("CHAMELEON_Lap2Desc Failed!");
    }
}

template<typename T>
void ChameleonImplementation<T>::ExaGeoStatLapackCopyTile(const common::UpperLower &aUpperLower, void *apA, void *apB) {
    int status = CHAMELEON_dlacpy_Tile((cham_uplo_t) aUpperLower, (CHAM_desc_t *) apA, (CHAM_desc_t *) apB);
    if (status != CHAMELEON_SUCCESS) {
        throw std::runtime_error("CHAMELEON_dlacpy_Tile Failed!");
    }
}

template<typename T>
void *ChameleonImplementation<T>::ExaGeoStatDataGetAddr(void *apA, int aAm, int aAn) {
    return RUNTIME_data_getaddr((CHAM_desc_t *) apA, aAm, aAn);
}

template<typename T>
void ChameleonImplementation<T>::ExaGeoStatTrsmTile(const common::Side &aSide, const common::UpperLower &aUpperLower,
                                                    const common::Trans &aTrans, const common::Diag &aDiag,
                                                    const T &aAlpha, void *apA, void *apCD, void *apCrk, void *apZ,
                                                    const int &aMaxRank) {
    int status = CHAMELEON_dtrsm_Tile((cham_side_t) aSide, (cham_uplo_t) aUpperLower, (cham_trans_t) aTrans,
                                      (cham_diag_t) aDiag, aAlpha, (CHAM_desc_t *) apA, (CHAM_desc_t *) apZ);
    if (status != CHAMELEON_SUCCESS) {
        throw std::runtime_error("CHAMELEON_dtrsm_Tile Failed!");
    }
}

template<typename T>
void ChameleonImplementation<T>::ExaGeoStatSequenceWait(void *apSequence) {
    int status = CHAMELEON_Sequence_Wait((RUNTIME_sequence_t *) apSequence);
    if (status != CHAMELEON_SUCCESS) {
        throw std::runtime_error("CHAMELEON_Sequence_Wait Failed!");
    }
}

template<typename T>
void ChameleonImplementation<T>::ExaGeoStatCreateSequence(void *apSequence) {
    int status = CHAMELEON_Sequence_Create(((RUNTIME_sequence_t **) apSequence));
    if (status != CHAMELEON_SUCCESS) {
        throw std::runtime_error("CHAMELEON_Sequence_Create Failed!");
    }
}

template<typename T>
void
ChameleonImplementation<T>::ExaGeoStatOptionsInit(void *apOptions, void *apContext, void *apSequence, void *apRequest) {
    RUNTIME_options_init((RUNTIME_option_t *) apOptions, (CHAM_context_t *) apContext,
                         (RUNTIME_sequence_t *) apSequence, (RUNTIME_request_t *) apRequest);
}

template<typename T>
void ChameleonImplementation<T>::ExaGeoStatOptionsFree(void *apOptions) {
    RUNTIME_options_ws_free((RUNTIME_option_t *) apOptions);
}

template<typename T>
void ChameleonImplementation<T>::ExaGeoStatOptionsFinalize(void *apOptions, void *apContext) {
    RUNTIME_options_finalize((RUNTIME_option_t *) apOptions, (CHAM_context_t *) apContext);
}

template<typename T>
int ChameleonImplementation<T>::ExaGeoStatMeasureDetTileAsync(void *apDescA, void *apSequence, void *apRequest,
                                                              void *apDescDet) {

    // Check for initialize the Chameleon context.
    if (!this->mpContext) {
        throw std::runtime_error(
                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStatHardware(computation, cores_number, gpu_numbers);'.");
    }
    RUNTIME_option_t options;
    this->ExaGeoStatOptionsInit(&options, this->mpContext, apSequence, apRequest);

    int m;
    int temp;
    auto Z = (CHAM_desc_t *) apDescA;
    auto det = (CHAM_desc_t *) apDescDet;
    struct starpu_codelet *cl = &this->cl_dmdet;

    for (m = 0; m < Z->mt; m++) {
        temp = m == Z->mt - 1 ? Z->m - m * Z->mb : Z->mb;
        starpu_insert_task(cl,
                           STARPU_VALUE, &temp, sizeof(int),
                           STARPU_R, ExaGeoStatDataGetAddr(Z, m, m),
                           STARPU_RW, ExaGeoStatDataGetAddr(det, 0, 0),
                           0);
    }
    this->ExaGeoStatOptionsFree(&options);
    this->ExaGeoStatOptionsFinalize(&options, (CHAM_context_t *) this->mpContext);
    return CHAMELEON_SUCCESS;
}

template<typename T>
int ChameleonImplementation<T>::ExaGeoStaStrideVectorTileAsync(void *apDescA, void *apDescB, void *apDescC,
                                                               void *apSequence, void *apRequest) {

    // Check for initialize the Chameleon context.
    if (!this->mpContext) {
        throw std::runtime_error(
                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStatHardware(computation, cores_number, gpu_numbers);'.");
    }

    RUNTIME_option_t options;
    this->ExaGeoStatOptionsInit(&options, this->mpContext, apSequence, apRequest);

    int m, m0;
    int temp;
    auto A = (CHAM_desc_t *) apDescA;
    auto B = (CHAM_desc_t *) apDescB;
    auto C = (CHAM_desc_t *) apDescC;
    struct starpu_codelet *cl = &this->cl_stride_vec;

    for (m = 0; m < A->mt; m++) {
        temp = m == A->mt - 1 ? A->m - m * A->mb : A->mb;
        m0 = m * A->mb;
        starpu_insert_task(cl,
                           STARPU_VALUE, &temp, sizeof(int),
                           STARPU_VALUE, &m0, sizeof(int),
                           STARPU_VALUE, &m, sizeof(int),
                           STARPU_R, ExaGeoStatDataGetAddr(A, m, 0),
                           STARPU_W, ExaGeoStatDataGetAddr(B, (int) floor(m / 2.0), 0),
                           STARPU_W, ExaGeoStatDataGetAddr(C, (int) floor(m / 2.0), 0),
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
                STARPU_NAME, "stride_vec",
#endif
                           0);

    }
    this->ExaGeoStatOptionsFree(&options);
    this->ExaGeoStatOptionsFinalize(&options, (CHAM_context_t *) this->mpContext);
    return CHAMELEON_SUCCESS;
}

template<typename T>
int
ChameleonImplementation<T>::ExaGeoStaStrideVectorTileAsync(void *apDescA, void *apDescB, void *apDescC, void *apDescD,
                                                           void *apSequence, void *apRequest) {

    // Check for initialize the Chameleon context.
    if (!this->mpContext) {
        throw std::runtime_error(
                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStatHardware(computation, cores_number, gpu_numbers);'.");
    }

    RUNTIME_option_t options;
    RUNTIME_options_init(&options, (CHAM_context_t *) this->mpContext, (RUNTIME_sequence_t *) apSequence,
                         (RUNTIME_request_t *) apRequest);

    int m, m0;
    int temp;
    auto A = (CHAM_desc_t *) apDescA;
    auto B = (CHAM_desc_t *) apDescB;
    auto C = (CHAM_desc_t *) apDescC;
    auto D = (CHAM_desc_t *) apDescD;
    struct starpu_codelet *cl = &this->cl_tri_stride_vec;

    for (m = 0; m < A->mt; m++) {
        temp = m == A->mt - 1 ? A->m - m * A->mb : A->mb;
        m0 = m * A->mb;
        starpu_insert_task(cl,
                           STARPU_VALUE, &temp, sizeof(int),
                           STARPU_VALUE, &m0, sizeof(int),
                           STARPU_VALUE, &m, sizeof(int),
                           STARPU_R, ExaGeoStatDataGetAddr(A, m, 0),
                           STARPU_W, ExaGeoStatDataGetAddr(B, (int) floor(m / 3.0), 0),
                           STARPU_W, ExaGeoStatDataGetAddr(C, (int) floor(m / 3.0), 0),
                           STARPU_W, ExaGeoStatDataGetAddr(D, (int) floor(m / 3.0), 0),
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
                STARPU_NAME, "tristride_vec",
#endif
                           0);

    }
    this->ExaGeoStatOptionsFree(&options);
    this->ExaGeoStatOptionsFinalize(&options, (CHAM_context_t *) this->mpContext);
    return CHAMELEON_SUCCESS;
}