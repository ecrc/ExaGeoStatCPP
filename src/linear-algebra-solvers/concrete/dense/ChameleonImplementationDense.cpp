
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ChameleonAllocateDescriptors.cpp
 * @brief Sets up the Chameleon descriptors needed for the dense matrix computations in ExaGeoStat.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-03-20
**/

#include <lapacke.h>

// Include Chameleon libraries
extern "C" {
#include <chameleon/struct.h>
#include <chameleon.h>
#include <control/descriptor.h>
#include <control/context.h>
}

#include <linear-algebra-solvers/concrete/dense/ChameleonImplementationDense.hpp>
#include <data-units/DescriptorData.hpp>

using namespace std;

using namespace exageostat::linearAlgebra::dense;
using namespace exageostat::common;
using namespace exageostat::dataunits;
using namespace exageostat::kernels;
using namespace exageostat::helpers;
using namespace exageostat::configurations;

// Define a method to set up the Chameleon descriptors
template<typename T>
void ChameleonImplementationDense<T>::InitiateDescriptors(Configurations *apConfigurations,
                                                          DescriptorData<T> *apDescriptorData) {

    // Check for Initialise the Chameleon context.
    if (!this->apContext) {
        throw std::runtime_error(
                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStat<double/float>::ExaGeoStatInitializeHardware(configurations)'.");
    }

    // Get the problem size and other configuration parameters
    int N = apConfigurations->GetProblemSize();
    int dts = apConfigurations->GetDenseTileSize();
    int p_grid = apConfigurations->GetPGrid();
    int q_grid = apConfigurations->GetQGrid();
    bool is_OOC = apConfigurations->GetIsOOC();

    // For distributed system and should be removed
    T *Zcpy = (T *) malloc(N * sizeof(T));
    T dot_product_value;

    // Create a Chameleon sequence
    RUNTIME_sequence_t *pSequence;
    RUNTIME_request_t request_array[2] = {CHAMELEON_SUCCESS, CHAMELEON_SUCCESS};
    CHAMELEON_Sequence_Create(&pSequence);


    // Set the floating point precision based on the template type
    FloatPoint float_point;
    if (sizeof(T) == SIZE_OF_FLOAT) {
        float_point = EXAGEOSTAT_REAL_FLOAT;
    } else if (sizeof(T) == SIZE_OF_DOUBLE) {
        float_point = EXAGEOSTAT_REAL_DOUBLE;
    } else {
        throw runtime_error("Unsupported for now!");
    }

//    double det = apConfigurations->GetDeterminantValue();

    apDescriptorData->SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_C, is_OOC, nullptr, float_point, dts, dts,
                                    dts * dts, N, N, 0, 0, N, N, p_grid, q_grid);
    apDescriptorData->SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_Z, is_OOC, nullptr, float_point, dts, dts,
                                    dts * dts, N, 1, 0, 0, N, 1, p_grid, q_grid);
    apDescriptorData->SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_Z_COPY, is_OOC, nullptr, float_point, dts,
                                    dts, dts * dts, N, 1, 0, 0, N, 1, p_grid, q_grid);
    apDescriptorData->SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_DETERMINANT, is_OOC, nullptr, float_point,
                                    dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);
    //// TODO: Modify this to be dot_product instead of det: nullptr works fine for poth det and prod, why?
    apDescriptorData->SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_PRODUCT, is_OOC, nullptr,
                                    float_point, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);

    if (float_point == EXAGEOSTAT_REAL_DOUBLE) {
        auto *CHAM_descC = apDescriptorData->GetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_C).chameleon_desc;
        apDescriptorData->CreateSubMatrixDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_C11, CHAM_descC, 0, 0,
                                                    CHAM_descC->m / 2, CHAM_descC->n / 2);
        apDescriptorData->CreateSubMatrixDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_C12, CHAM_descC,
                                                    CHAM_descC->m / 2, 0, CHAM_descC->m / 2, CHAM_descC->n / 2);
        apDescriptorData->CreateSubMatrixDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_C22, CHAM_descC,
                                                    CHAM_descC->m / 2, CHAM_descC->n / 2, CHAM_descC->m / 2,
                                                    CHAM_descC->n / 2);

        apDescriptorData->SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_Z_1, is_OOC, nullptr, float_point, dts,
                                        dts, dts * dts, N / 2, 1, 0, 0, N / 2, 1, p_grid, q_grid);
        apDescriptorData->SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_Z_2, is_OOC, nullptr, float_point, dts,
                                        dts, dts * dts, N / 2, 1, 0, 0, N / 2, 1, p_grid, q_grid);
        apDescriptorData->SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_PRODUCT_1, is_OOC, &dot_product_value,
                                        float_point, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);
        apDescriptorData->SetDescriptor(common::CHAMELEON_DESCRIPTOR, DESCRIPTOR_PRODUCT_2, is_OOC, &dot_product_value,
                                        float_point, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);
    }

    apDescriptorData->SetSequence(pSequence);
    apDescriptorData->SetRequest(request_array);

    //stop gsl error handler
    gsl_set_error_handler_off();
    free(Zcpy);
}

template<typename T>
void ChameleonImplementationDense<T>::ExaGeoStatInitContext(int aCoresNumber, int aGPUsNumbers) {

    if (!this->apContext) {
        CHAMELEON_user_tag_size(31, 26);
        CHAMELEON_Init(aCoresNumber, aGPUsNumbers)
        this->apContext = chameleon_context_self();
    }
}

template<typename T>
void ChameleonImplementationDense<T>::ExaGeoStatFinalizeContext() {

    if (!this->apContext) {
        cout
                << "No initialised context of Chameleon, Please use 'ExaGeoStat<double/or/float>::ExaGeoStatInitializeHardware(configurations);'"
                << endl;
    } else {
        CHAMELEON_Finalize()
        this->apContext = nullptr;
    }
}

template<typename T>
void ChameleonImplementationDense<T>::CovarianceMatrixCodelet(DescriptorData<T> *apDescriptorData, void *apDescriptor,
                                                              int &aTriangularPart,
                                                              dataunits::Locations<T> *apLocation1,
                                                              dataunits::Locations<T> *apLocation2,
                                                              dataunits::Locations<T> *apLocation3,
                                                              T *aLocalTheta, int aDistanceMetric,
                                                              kernels::Kernel<T> *apKernel) {

    // Check for Initialise the Chameleon context.
    if (!this->apContext) {
        throw std::runtime_error(
                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStat<double/float>::ExaGeoStatInitializeHardware(configurations)'.");
    }

    RUNTIME_option_t options;
    RUNTIME_options_init(&options, (CHAM_context_t *) this->apContext,
                         (RUNTIME_sequence_t *) apDescriptorData->GetSequence(),
                         (RUNTIME_request_t *) apDescriptorData->GetRequest());

    int tempmm, tempnn;
    auto *CHAM_apDescriptor = (CHAM_desc_t *) apDescriptor;
    CHAM_desc_t A = *CHAM_apDescriptor;

    struct starpu_codelet *cl = &this->cl_dcmg;
    int m, n, m0 = 0, n0 = 0;
    for (n = 0; n < A.nt; n++) {
        tempnn = n == A.nt - 1 ? A.n - n * A.nb : A.nb;
        if (aTriangularPart == ChamUpperLower) {
            m = 0;
        } else {
            m = A.m == A.n ? n : 0;
        }
        for (; m < A.mt; m++) {

            tempmm = m == A.mt - 1 ? A.m - m * A.mb : A.mb;
            m0 = m * A.mb;
            n0 = n * A.nb;

            // Register the data with StarPU
            starpu_insert_task(cl,
                               STARPU_VALUE, &tempmm, sizeof(int),
                               STARPU_VALUE, &tempnn, sizeof(int),
                               STARPU_VALUE, &m0, sizeof(int),
                               STARPU_VALUE, &n0, sizeof(int),
                               STARPU_W, (starpu_data_handle_t) RUNTIME_data_getaddr(CHAM_apDescriptor, m, n),
                               STARPU_VALUE, &apLocation1, sizeof(dataunits::Locations<T> *),
                               STARPU_VALUE, &apLocation2, sizeof(dataunits::Locations<T> *),
                               STARPU_VALUE, &apLocation3, sizeof(dataunits::Locations<T> *),
                               STARPU_VALUE, &aLocalTheta, sizeof(double *),
                               STARPU_VALUE, &aDistanceMetric, sizeof(int),
                               STARPU_VALUE, &apKernel, sizeof(kernels::Kernel<T> *),
                               0);
        }
    }
    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, (CHAM_context_t *) this->apContext);

    CHAMELEON_Sequence_Wait((RUNTIME_sequence_t *) apDescriptorData->GetSequence());
}

template<typename T>
void ChameleonImplementationDense<T>::GenerateObservationsVector(Configurations *apConfigurations,
                                                                 DescriptorData<T> *apDescriptorData,
                                                                 BaseDescriptor aDescriptor, Locations<T> *apLocation1,
                                                                 Locations<T> *apLocation2, Locations<T> *apLocation3,
                                                                 int aDistanceMetric, Kernel<T> *apKernel, T *Nrand) {

    // Check for Initialise the Chameleon context.
    if (!this->apContext) {
        throw std::runtime_error(
                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStat<double/float>::ExaGeoStatInitializeHardware(configurations)'.");
    }

    const int N = apConfigurations->GetProblemSize();
    int seed = apConfigurations->GetSeed();
    int iseed[4] = {seed, seed, seed, 1};
    auto *pDescriptor = aDescriptor.chameleon_desc;

    //Generate the co-variance matrix C
    auto *theta = (T *) malloc(apConfigurations->GetInitialTheta().size() * sizeof(T));
    for (int i = 0; i < apConfigurations->GetInitialTheta().size(); i++) {
        theta[i] = apConfigurations->GetInitialTheta()[i];
    }

    VERBOSE("Initializing Covariance Matrix (Synthetic Dataset Generation Phase).....")
    int upper_lower = EXAGEOSTAT_LOWER;
    this->CovarianceMatrixCodelet(apDescriptorData, pDescriptor, upper_lower, apLocation1, apLocation2, apLocation3,
                                  theta,
                                  aDistanceMetric, apKernel);
    free(theta);
    VERBOSE("Done.\n")

    //Copy Nrand to Z
    VERBOSE("Generate Normal Random Distribution Vector Z (Synthetic Dataset Generation Phase) .....")
    auto *CHAM_descriptorZ = apDescriptorData->GetDescriptor(CHAMELEON_DESCRIPTOR, DESCRIPTOR_Z).chameleon_desc;
    CopyDescriptorZ(apDescriptorData, (void *) CHAM_descriptorZ, Nrand);
    VERBOSE("Done.\n")

    //Cholesky factorization for the Co-variance matrix C
    VERBOSE("Cholesky factorization of Sigma (Synthetic Dataset Generation Phase) .....")

    int potential_failure = CHAMELEON_dpotrf_Tile(ChamLower, pDescriptor);
    FAILURE_LOGGER(potential_failure, "Factorization cannot be performed..\nThe matrix is not positive definite")
    VERBOSE("Done.\n")

    //Triangular matrix-matrix multiplication
    VERBOSE("Triangular matrix-matrix multiplication Z=L.e (Synthetic Dataset Generation Phase) .....")
    CHAMELEON_dtrmm_Tile(ChamLeft, ChamLower, ChamNoTrans, ChamNonUnit, 1, pDescriptor, CHAM_descriptorZ);
    VERBOSE("Done.\n")

    const int P = apConfigurations->GetP();
    if (apConfigurations->GetLogger()) {
        T *pMatrix;
        VERBOSE("Writing generated data to the disk (Synthetic Dataset Generation Phase) .....")
#ifdef CHAMELEON_USE_MPI
        pMatrix = (T*) malloc(N * sizeof(T));
        CHAMELEON_Tile_to_Lapack( CHAM_descriptorZ, pMatrix, N);
        if ( CHAMELEON_My_Mpi_Rank() == 0 ){
            DiskWriter<T>::WriteVectorsToDisk(pMatrix, &N, &P, apConfigurations->GetLoggerPath(), apLocation1);
        }
        free(pMatrix);
#else
        pMatrix = (T *) CHAM_descriptorZ->mat;
        DiskWriter<T>::WriteVectorsToDisk(pMatrix, &N, &P, apConfigurations->GetLoggerPath(), apLocation1);
        free(pMatrix);
#endif
        VERBOSE(" Done.\n")
    }

    CHAMELEON_dlaset_Tile(ChamUpperLower, 0, 0, pDescriptor);
    free(Nrand);
    VERBOSE("Done Z Vector Generation Phase. (Chameleon Synchronous)")
}

template<typename T>
void
ChameleonImplementationDense<T>::CopyDescriptorZ(DescriptorData<T> *apDescriptorData, void *apDescriptor,
                                                 T *apDoubleVector) {

    // Check for Initialise the Chameleon context.
    if (!this->apContext) {
        throw std::runtime_error(
                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStat<double/float>::ExaGeoStatInitializeHardware(configurations)'.");
    }

    RUNTIME_option_t options;
    RUNTIME_options_init(&options, (CHAM_context_t *) this->apContext,
                         (RUNTIME_sequence_t *) apDescriptorData->GetSequence(),
                         (RUNTIME_request_t *) apDescriptorData->GetRequest());

    int m, m0;
    int tempmm;
    auto A = (CHAM_desc_t *) apDescriptor;
    struct starpu_codelet *cl = &this->cl_dzcpy;

    for (m = 0; m < A->mt; m++) {
        tempmm = m == A->mt - 1 ? A->m - m * A->mb : A->mb;
        m0 = m * A->mb;

        starpu_insert_task(cl,
                           STARPU_VALUE, &tempmm, sizeof(int),
                           STARPU_VALUE, &m0, sizeof(int),
                           STARPU_VALUE, &apDoubleVector, sizeof(double),
                           STARPU_W, RUNTIME_data_getaddr(A, m, 0),
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
                STARPU_NAME, "dzcpy",
#endif
                           0);
    }
    RUNTIME_options_ws_free(&options);

}

template<typename T>
void ChameleonImplementationDense<T>::DestroyDescriptors(DescriptorData<T> *apDescriptorData) {

    CHAM_desc_t *pDescriptorC = apDescriptorData->GetDescriptor(CHAMELEON_DESCRIPTOR, DESCRIPTOR_C).chameleon_desc;
    CHAM_desc_t *pDescriptorZ = apDescriptorData->GetDescriptor(CHAMELEON_DESCRIPTOR, DESCRIPTOR_Z).chameleon_desc;
    CHAM_desc_t *pDescriptorZ_1 = apDescriptorData->GetDescriptor(CHAMELEON_DESCRIPTOR, DESCRIPTOR_Z_1).chameleon_desc;
    CHAM_desc_t *pDescriptorZ_2 = apDescriptorData->GetDescriptor(CHAMELEON_DESCRIPTOR, DESCRIPTOR_Z_2).chameleon_desc;
    CHAM_desc_t *pDescriptorZ_COPY = apDescriptorData->GetDescriptor(CHAMELEON_DESCRIPTOR,
                                                                     DESCRIPTOR_Z_COPY).chameleon_desc;
    CHAM_desc_t *pDescriptorProduct = apDescriptorData->GetDescriptor(CHAMELEON_DESCRIPTOR,
                                                                      DESCRIPTOR_PRODUCT).chameleon_desc;
    CHAM_desc_t *pDescriptorProduct_1 = apDescriptorData->GetDescriptor(CHAMELEON_DESCRIPTOR,
                                                                        DESCRIPTOR_PRODUCT_1).chameleon_desc;
    CHAM_desc_t *pDescriptorProduct_2 = apDescriptorData->GetDescriptor(CHAMELEON_DESCRIPTOR,
                                                                        DESCRIPTOR_PRODUCT_2).chameleon_desc;
    CHAM_desc_t *pDescriptorDeterminant = apDescriptorData->GetDescriptor(CHAMELEON_DESCRIPTOR,
                                                                          DESCRIPTOR_DETERMINANT).chameleon_desc;

    if (pDescriptorC) {
        CHAMELEON_Desc_Destroy(&pDescriptorC);
    }
    if (pDescriptorZ) {
        CHAMELEON_Desc_Destroy(&pDescriptorZ);
    }
    if (pDescriptorZ_1) {
        CHAMELEON_Desc_Destroy(&pDescriptorZ_1);
    }
    if (pDescriptorZ_2) {
        CHAMELEON_Desc_Destroy(&pDescriptorZ_2);
    }
    if (pDescriptorZ_COPY) {
        CHAMELEON_Desc_Destroy(&pDescriptorZ_COPY);
    }
    if (pDescriptorProduct) {
        CHAMELEON_Desc_Destroy(&pDescriptorProduct);
    }
    if (pDescriptorProduct_1) {
        CHAMELEON_Desc_Destroy(&pDescriptorProduct_1);
    }
    if (pDescriptorProduct_2) {
        CHAMELEON_Desc_Destroy(&pDescriptorProduct_2);
    }
    if (pDescriptorDeterminant) {
        CHAMELEON_Desc_Destroy(&pDescriptorDeterminant);
    }

    if ((RUNTIME_sequence_t *) apDescriptorData->GetSequence()) {
        CHAMELEON_Sequence_Destroy((RUNTIME_sequence_t *) apDescriptorData->GetSequence());
    }
}

template<typename T>
void ChameleonImplementationDense<T>::ExaGeoStatAllocateMatrixTile(void **apDescriptor, bool ais_OOC, T *apMemSpace,
                                                                   int aType2, int aMB,
                                                                   int aNB, int aMBxNB, int aLda, int aN, int aSMB,
                                                                   int aSNB, int aM, int aN2, int aP, int aQ) {
    if (ais_OOC && apMemSpace == nullptr && aMB != 1 && aNB != 1) {
        CHAMELEON_Desc_Create_OOC((CHAM_desc_t **) apDescriptor, (cham_flttype_t) aType2, aMB, aNB, aMBxNB, aLda, aN,
                                  aSMB, aSNB, aM, aN2, aP, aQ);
    } else {
        CHAMELEON_Desc_Create((CHAM_desc_t **) apDescriptor, apMemSpace, (cham_flttype_t) aType2, aMB, aNB, aMBxNB,
                              aLda, aN, aSMB, aSNB, aM, aN2, aP, aQ);
    }
}

template<typename T>
T ChameleonImplementationDense<T>::ExaGeoStatMleTile(ExaGeoStatData<T> *apData,
                                                     configurations::Configurations *apConfigurations,
                                                     const double *theta) {
    cout << "\t\t **MLE**\n";

    //Initialization
    T loglik = 0.0, logdet = 0.0, time_facto = 0.0, time_solve = 0.0, logdet_calculate = 0.0, matrix_gen_time = 0.0;
    T dzcpy_time = 0.0, variance = 1, variance1 = 1, variance2 = 1, variance3 = 1, dot_product = 0, dot_product1 = 0,
    dot_product2 = 0, dot_product3 = 0;

    double avg_executed_time_per_iteration = 0, avg_flops_per_iter = 0.0;

    int n, nhrs, success, i, num_params;
    T flops = 0.0;
    T *univariate_theta;
    T *univariate2_theta;
    T *univariate3_theta;
    T nu12;
    T rho;
    T sigma_square12;

    auto kernel_name = apConfigurations->GetKernelName();
    auto median_locations = apData->CalculateMedianLocations(kernel_name);

    auto *CHAM_desc_C = apData->GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                   DescriptorName::DESCRIPTOR_C).chameleon_desc;
    auto *CHAM_desc_sub_C11 = apData->GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                         DescriptorName::DESCRIPTOR_C11).chameleon_desc;
    auto *CHAM_desc_sub_C12 = apData->GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                         DescriptorName::DESCRIPTOR_C12).chameleon_desc;
    auto *CHAM_desc_sub_C22 = apData->GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                         DescriptorName::DESCRIPTOR_C22).chameleon_desc;

    auto *CHAM_desc_Z = apData->GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                   DescriptorName::DESCRIPTOR_Z).chameleon_desc;
    auto *CHAM_desc_Z1 = apData->GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                    DescriptorName::DESCRIPTOR_Z_1).chameleon_desc;
    auto *CHAM_desc_Z2 = apData->GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                    DescriptorName::DESCRIPTOR_Z_2).chameleon_desc;
    auto *CHAM_desc_Zcpy = apData->GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                      DescriptorName::DESCRIPTOR_Z_COPY).chameleon_desc;

    auto *CHAM_desc_det = apData->GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                     DescriptorName::DESCRIPTOR_DETERMINANT).chameleon_desc;

    auto *CHAM_desc_product = apData->GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                         DescriptorName::DESCRIPTOR_PRODUCT).chameleon_desc;
    auto *CHAM_desc_product1 = apData->GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                          DescriptorName::DESCRIPTOR_PRODUCT_1).chameleon_desc;
    auto *CHAM_desc_product2 = apData->GetDescriptorData()->GetDescriptor(DescriptorType::CHAMELEON_DESCRIPTOR,
                                                                          DescriptorName::DESCRIPTOR_PRODUCT_2).chameleon_desc;

    auto sequence = (RUNTIME_sequence_t *) apData->GetDescriptorData()->GetSequence();
    auto request = (RUNTIME_request_t *) apData->GetDescriptorData()->GetRequest();

    T *determinant = apData->GetDescriptorData()->GetDescriptorMatrix(CHAM_desc_det);
    *determinant = 0;

    T *product = apData->GetDescriptorData()->GetDescriptorMatrix(CHAM_desc_product);
    *product = 0;

    n = CHAM_desc_C->m;
    nhrs = CHAM_desc_Z->n;

    START_TIMING(dzcpy_time);
    int iter_count = apConfigurations->GetIterationsValue();
    if (iter_count == 0) {
    //Save a copy of descZ into descZcpy for restoring each iteration (Only for the first iteration)
        ExaGeoStatLapackCopyTile(UpperLower::EXAGEOSTAT_UPPER_LOWER, CHAM_desc_Z, CHAM_desc_Zcpy);
    }
    string recovery_file = apConfigurations->GetRecoveryFile();

    if (recovery_file.empty() ||
        !(this->recover((char *) (recovery_file.c_str()), iter_count, (T *) theta, &loglik, num_params))) {
        START_TIMING(dzcpy_time);
        if (iter_count == 0) {
    //    Save a copy of descZ into descZcpy for restoring each iteration (Only for the first iteration)
        ExaGeoStatLapackCopyTile(EXAGEOSTAT_UPPER_LOWER, CHAM_desc_Z, CHAM_desc_Zcpy);
        } else {
            VERBOSE("Re-store the original Z vector...");
            ExaGeoStatLapackCopyTile(EXAGEOSTAT_UPPER_LOWER, CHAM_desc_Zcpy, CHAM_desc_Z);
            VERBOSE(" Done.\n")
        }
        STOP_TIMING(dzcpy_time)

        //Generate new co-variance matrix C based on new theta
        VERBOSE("Generate New Covariance Matrix...")
        START_TIMING(matrix_gen_time)

        if (kernel_name == "bivariate_matern_parsimonious2" ||
            kernel_name == "bivariate_matern_parsimonious2_profile") {

            int upper_lower = EXAGEOSTAT_UPPER_LOWER;
            kernels::Kernel<T> *kernel = exageostat::plugins::PluginRegistry<kernels::Kernel<T>>::Create(
                    kernel_name);
            num_params = kernel->GetParametersNumbers();
            univariate_theta = (T *) malloc(3 * sizeof(T));
            univariate2_theta = (T *) malloc(3 * sizeof(T));
            univariate3_theta = (T *) malloc(3 * sizeof(T));
            univariate_theta[0] = theta[0];
            univariate_theta[1] = theta[2];
            univariate_theta[2] = theta[3];

            int distance_metric = apConfigurations->GetDistanceMetric();

            this->CovarianceMatrixCodelet(apData->GetDescriptorData(), CHAM_desc_sub_C11, upper_lower,
                                          apData->GetLocations(), apData->GetLocations(),
                                          median_locations, univariate_theta,
                                          distance_metric, kernel);

            nu12 = 0.5 * (theta[3] + theta[4]);

            rho = theta[5] * sqrt((tgamma(theta[3] + 1) * tgamma(theta[4] + 1)) /
                                  (tgamma(theta[3]) * tgamma(theta[4]))) *
                  tgamma(nu12) / tgamma(nu12 + 1);
            sigma_square12 = rho * sqrt(theta[0] * theta[1]);

            univariate2_theta[0] = sigma_square12;
            univariate2_theta[1] = theta[2];
            univariate2_theta[2] = nu12;

            this->CovarianceMatrixCodelet(apData->GetDescriptorData(), CHAM_desc_sub_C12, upper_lower, median_locations,
                                          apData->GetLocations(),
                                          median_locations, univariate2_theta,
                                          0, kernel);

            STOP_TIMING(matrix_gen_time);
            VERBOSE(" Done.\n");

            univariate3_theta[0] = theta[1];
            univariate3_theta[1] = theta[2];
            univariate3_theta[2] = theta[4];

            this->CovarianceMatrixCodelet(apData->GetDescriptorData(), CHAM_desc_sub_C22, upper_lower, median_locations,
                                          apData->GetLocations(),
                                          median_locations, univariate2_theta,
                                          0, kernel);
    } else {
        int upper_lower = EXAGEOSTAT_LOWER;

        kernels::Kernel<T> *kernel = exageostat::plugins::PluginRegistry<kernels::Kernel<T>>::Create(
                kernel_name);
        num_params = kernel->GetParametersNumbers();

        ((T *) apData->GetLocations()->GetLocationY())[0] = 0.103883;
        ((T *) apData->GetLocations()->GetLocationY())[1] = 0.135790;

        this->CovarianceMatrixCodelet(apData->GetDescriptorData(), CHAM_desc_C, upper_lower, apData->GetLocations(),
                                      apData->GetLocations(),
                                      apData->GetLocations(), (T *) theta,
                                      0, kernel);
    }

    ExaGeoStatSequenceWait(sequence);
    STOP_TIMING(matrix_gen_time);

    VERBOSE("Cholesky factorization of Sigma...");
    START_TIMING(time_facto);

    success = CHAMELEON_dpotrf_Tile(ChamLower, CHAM_desc_C);
    FAILURE_LOGGER(success, "Factorization cannot be performed..\n The matrix is not positive definite\n\n");
    STOP_TIMING(time_facto);
    flops = flops + FLOPS_DPOTRF(n);
    VERBOSE(" Done.\n");

    //Calculate log(|C|) --> log(square(|L|))
    VERBOSE("Calculating the log determinant ...");
    START_TIMING(logdet_calculate);
    ExaGeoStatMeasureDetTileAsync(CHAM_desc_C, sequence, &request, CHAM_desc_det);
    ExaGeoStatSequenceWait(sequence);

    //TODO: determinant value get rounded up, which affects results.
    logdet = 2 * (*determinant);
    STOP_TIMING(logdet_calculate);
    VERBOSE(" Done.\n");

//    Solving Linear System (L*X=Z)--->inv(L)*Z
    VERBOSE("Solving the linear system ...\n");
    START_TIMING(time_solve);
    ExaGeoStatTrsmTile(EXAGEOSTAT_LEFT, EXAGEOSTAT_LOWER, EXAGEOSTAT_NO_TRANS, EXAGEOSTAT_NON_UNIT, 1, CHAM_desc_C,
                       CHAM_desc_Z);
    STOP_TIMING(time_solve);
    flops = flops + FLOPS_DTRSM(ChamLeft, n, nhrs);
    VERBOSE(" Done.\n")

    //Calculate MLE likelihood
    VERBOSE("Calculating the MLE likelihood function ...");

    ExaGeoStatGemmTile(EXAGEOSTAT_TRANS, EXAGEOSTAT_NO_TRANS, 1, CHAM_desc_Z, CHAM_desc_Z, 0, CHAM_desc_product);
    if (kernel_name == "bivariate_matern_parsimonious2_profile") {

        loglik = -(n / 2) + (n / 2) * log(n) - (n / 2) * log(dot_product) - 0.5 * logdet -
                 (T) (n / 2.0) * log(2.0 * PI);
        ExaGeoStatGemmTile(EXAGEOSTAT_TRANS, EXAGEOSTAT_NO_TRANS, 1, CHAM_desc_Z1, CHAM_desc_Z1, 0,
                           CHAM_desc_product1);
        ExaGeoStatGemmTile(EXAGEOSTAT_TRANS, EXAGEOSTAT_NO_TRANS, 1, CHAM_desc_Z2, CHAM_desc_Z2, 0,
                           CHAM_desc_product2);

        variance1 = (1.0 / (n / 2)) * dot_product1;
        variance2 = (1.0 / (n / 2)) * dot_product2;

    } else if (kernel_name == "bivariate_matern_parsimonious_profile") {
        loglik = -(n / 2) + (n / 2) * log(n) - (n / 2) * log(dot_product) - 0.5 * logdet -
                 (double) (n / 2.0) * log(2.0 * PI);

        //TODO: STRIDE

        ExaGeoStatGemmTile(EXAGEOSTAT_TRANS, EXAGEOSTAT_NO_TRANS, 1, CHAM_desc_Z1, CHAM_desc_Z1, 0,
                           CHAM_desc_product1);
        ExaGeoStatGemmTile(EXAGEOSTAT_TRANS, EXAGEOSTAT_NO_TRANS, 1, CHAM_desc_Z2, CHAM_desc_Z2, 0,
                           CHAM_desc_product2);
        variance1 = (1.0 / (n / 2)) * dot_product1;
        variance2 = (1.0 / (n / 2)) * dot_product2;

    } else {
        dot_product = *product;
        loglik = -0.5 * dot_product - 0.5 * logdet - (double) (n / 2.0) * log(2.0 * PI);
        variance = theta[0];
    }
      VERBOSE(" Done.\n");
    }

    fprintf(stderr, " %3d- Model Parameters (", iter_count + 1);

    if (apConfigurations->GetLogger()) {
        ///should we be leaving the file opening and closing to the user?
        fprintf(apConfigurations->GetFileLog(), " %3d- Model Parameters (", iter_count + 1);
    }

    if ((apConfigurations->GetKernelName() == "bivariate_matern_parsimonious_profile") ||
        (apConfigurations->GetKernelName() == "bivariate_matern_parsimonious2_profile")) {
        fprintf(stderr, "%.8f, %.8f,", variance1, variance2);
        if (apConfigurations->GetLogger()) {
            fprintf(apConfigurations->GetFileLog(), "%.8f, %.8f,", variance1,
                    variance2);
        }
        i = 2;
    } else {
        i = 0;
    }

    for (; i < num_params; i++) {
        fprintf(stderr, "%.8f", theta[i]);
        if (i < num_params - 1) {
            fprintf(stderr, ",");
        }

        if (apConfigurations->GetLogger()) {
            fprintf(apConfigurations->GetFileLog(), "%.8f, ", theta[i]);
        }
    }
    ExaGeoStatSequenceWait(sequence);
    fprintf(stderr, ")----> LogLi: %.18f\n", loglik);

    if (apConfigurations->GetLogger()) {
        fprintf(apConfigurations->GetFileLog(), ")----> LogLi: %.18f\n", loglik);
    }


    fprintf(stderr, " ---- Facto Time: %6.2f \n", time_facto);
    fprintf(stderr, " ---- Matrix Generation Time: %6.2f\n", matrix_gen_time);
    fprintf(stderr, " ---- Total Time: %6.2f\n", matrix_gen_time + time_facto + logdet_calculate + time_solve);

    apConfigurations->SetIterationsValue(apConfigurations->GetIterationsValue() + 1);
    // for experiments
    if (Configurations::GetRunMode() == VERBOSE_MODE) {
        avg_executed_time_per_iteration += /*matrix_gen_time*/+time_facto + logdet_calculate + time_solve;
        avg_flops_per_iter += flops / 1e9 / (time_facto + time_solve);
    }
    apConfigurations->SetFinalLogLik(loglik);

//    output
    return loglik;
}

template<typename T>
int ChameleonImplementationDense<T>::ExaGeoStatLapackCopyTile(exageostat::common::UpperLower aUpperLower, void *apA,
                                                              void *apB) {
    return CHAMELEON_dlacpy_Tile((cham_uplo_t) aUpperLower, (CHAM_desc_t *) apA, (CHAM_desc_t *) apB);
}

template<typename T>
int
ChameleonImplementationDense<T>::ExaGeoStatLapackToDescriptor(exageostat::common::UpperLower aUpperLower, void *apAf77,
                                                              int aLda, void *apA) {
    return CHAMELEON_Lap2Desc((cham_uplo_t) aUpperLower, (CHAM_desc_t *) apAf77, aLda, (CHAM_desc_t *) apA);
}


template<typename T>
int ChameleonImplementationDense<T>::ExaGeoStatSequenceWait(void *apSequence) {
    return CHAMELEON_Sequence_Wait((RUNTIME_sequence_t *) apSequence);
}

template<typename T>
int ChameleonImplementationDense<T>::ExaGeoStatPotrfTile(exageostat::common::UpperLower aUpperLower, void *apA) {
    return CHAMELEON_dpotrf_Tile((cham_uplo_t) aUpperLower, (CHAM_desc_t *) apA);
}

template<typename T>
int ChameleonImplementationDense<T>::ExaGeoStatTrsmTile(common::Side aSide, common::UpperLower aUpperLower,
                                                        common::Trans aTrans, common::Diag aDiag, T aAlpha, void *apA,
                                                        void *apB) {
    return CHAMELEON_dtrsm_Tile((cham_side_t) aSide, (cham_uplo_t) aUpperLower, (cham_trans_t) aTrans,
                                (cham_diag_t) aDiag, aAlpha, (CHAM_desc_t *) apA, (CHAM_desc_t *) apB);
}

template<typename T>
int ChameleonImplementationDense<T>::ExaGeoStatGemmTile(common::Trans aTransA, common::Trans aTransB, T aAlpha,
                                                        void *apA, void *apB, T aBeta, void *apC) {
    return CHAMELEON_dgemm_Tile((cham_trans_t) aTransA, (cham_trans_t) aTransB, aAlpha, (CHAM_desc_t *) apA,
                                (CHAM_desc_t *) apB, aBeta, (CHAM_desc_t *) apC);
}

template<typename T>
int ChameleonImplementationDense<T>::ExaGeoStatMeasureDetTileAsync(void *apDescA, void *apSequence, void *apRequest,
                                                                   void *apDescDet) {

    // Check for Initialise the Chameleon context.
    if (!this->apContext) {
        throw std::runtime_error(
                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStat<double/float>::ExaGeoStatInitializeHardware(configurations)'.");
    }
    RUNTIME_option_t options;
    RUNTIME_options_init(&options, (CHAM_context_t *) this->apContext,
                         (RUNTIME_sequence_t *) apSequence, (RUNTIME_request_t *) apRequest);

    int m;
    int tempmm;
    auto Z = (CHAM_desc_t *) apDescA;
    CHAM_desc_t A = *Z;
    auto det = (CHAM_desc_t *) apDescDet;
    struct starpu_codelet *cl = &this->cl_dmdet;

    for (m = 0; m < A.mt; m++) {
        tempmm = m == A.mt - 1 ? A.m - m * A.mb : A.mb;

        starpu_insert_task(cl,
                           STARPU_VALUE, &tempmm, sizeof(int),
                           STARPU_R, RUNTIME_data_getaddr(Z, m, m),
                           STARPU_RW, RUNTIME_data_getaddr(det, 0, 0),
                           0);
    }
    RUNTIME_options_ws_free(&options);
    return CHAMELEON_SUCCESS;
}


namespace exageostat::linearAlgebra::dense {
    template<typename T> void *ChameleonImplementationDense<T>::apContext = nullptr;
}


template<typename T>
int ChameleonImplementationDense<T>::ExaGeoStaStrideVectorTileAsync(void *apDescA, void *apDescB, void *apDescC,
                                                                    void *apSequence, void *apRequest) {
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    chamctxt = chameleon_context_self();
    if (((RUNTIME_sequence_t *) apSequence)->status != CHAMELEON_SUCCESS)
        return -2;
    RUNTIME_options_init(&options, chamctxt, (RUNTIME_sequence_t *) apSequence, (RUNTIME_request_t *) apRequest);

    int m, m0;
    int tempmm;
    CHAM_desc_t A = *((CHAM_desc_t *) apDescA);
    CHAM_desc_t B = *((CHAM_desc_t *) apDescB);
    CHAM_desc_t C = *((CHAM_desc_t *) apDescC);
    struct starpu_codelet *cl = &this->cl_stride_vec;
    for (m = 0; m < A.mt; m++) {
        tempmm = m == A.mt - 1 ? A.m - m * A.mb : A.mb;
        m0 = m * A.mb;
        starpu_insert_task(cl,
                           STARPU_VALUE, &tempmm, sizeof(int),
                           STARPU_VALUE, &m0, sizeof(int),
                           STARPU_VALUE, &m, sizeof(int),
                           STARPU_R, (starpu_data_handle_t) RUNTIME_data_getaddr((CHAM_desc_t *) apDescA, m, 0),
                           STARPU_W,
                           (starpu_data_handle_t) RUNTIME_data_getaddr((CHAM_desc_t *) apDescB, (int) floor(m / 2.0),
                                                                       0),
                           STARPU_W,
                           (starpu_data_handle_t) RUNTIME_data_getaddr((CHAM_desc_t *) apDescC, (int) floor(m / 2.0),
                                                                       0),
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
                STARPU_NAME, "stride_vec",
#endif
                           0);

    }
    RUNTIME_options_ws_free(&options);
    return CHAMELEON_SUCCESS;
}
