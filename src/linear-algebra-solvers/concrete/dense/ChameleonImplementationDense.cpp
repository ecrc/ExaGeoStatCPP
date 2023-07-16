
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
#include <configurations/data-generation/concrete/SyntheticDataConfigurations.hpp>

using namespace std;

using namespace exageostat::linearAlgebra::dense;
using namespace exageostat::common;
using namespace exageostat::dataunits;
using namespace exageostat::kernels;
using namespace exageostat::helpers;
using namespace exageostat::configurations::data_modeling;

// Define a method to set up the Chameleon descriptors
template<typename T>
void ChameleonImplementationDense<T>::InitiateDescriptors() {

    // Check for Initialise the Chameleon context.
    if (!this->apContext) {
        throw std::runtime_error(
                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStat<double/float>::ExaGeoStatInitializeHardware(configurations)'.");
    }

    // Declare variables for Chameleon descriptors
    vector<void *> &pDescriptorC = this->mpConfigurations->GetDescriptorC();
    vector<void *> &pDescriptorZ = this->mpConfigurations->GetDescriptorZ();
    auto pDescriptorZcpy = &this->mpConfigurations->GetDescriptorZcpy();
    vector<void *> &pDescriptorProduct = this->mpConfigurations->GetDescriptorProduct();
    auto pDescriptorDeterminant = &this->mpConfigurations->GetDescriptorDeterminant();

    // Get the problem size and other configuration parameters
    int vector_size;
    int N = this->mpConfigurations->GetProblemSize();
    int dts = this->mpConfigurations->GetDenseTileSize();
    int p_grid = this->mpConfigurations->GetPGrid();
    int q_grid = this->mpConfigurations->GetQGrid();
    bool is_OOC = this->mpConfigurations->GetIsOOC();

    // For distributed system and should be removed
    T *Zcpy = (T *) malloc(N * sizeof(T));
    T dot_product_value;

    // Create a Chameleon sequence
    RUNTIME_sequence_t *pSequence;
    RUNTIME_request_t request[2] = {CHAMELEON_SUCCESS, CHAMELEON_SUCCESS};
    CHAMELEON_Sequence_Create(&pSequence);

    // Set the floating point precision based on the template type
    FloatPoint float_point;
    if (sizeof(T) == SIZE_OF_FLOAT) {
        float_point = EXAGEOSTAT_REAL_FLOAT;
        vector_size = 1;
    } else {
        float_point = EXAGEOSTAT_REAL_DOUBLE;
        vector_size = 3;
    }

    // Create the Chameleon descriptors based on the configuration parameters
    // Depending on the passed Precession, the descriptor will resize for value 1 or 3.
    pDescriptorC.resize(vector_size + 1, nullptr);
    pDescriptorZ.resize(vector_size, nullptr);
    pDescriptorProduct.resize(vector_size, nullptr);

    auto **CHAM_descriptorC = &pDescriptorC[0];
    ExageostatAllocateMatrixTile(CHAM_descriptorC, is_OOC, nullptr, (cham_flttype_t) float_point, dts, dts,
                                 dts * dts, N, N, 0, 0, N, N, p_grid, q_grid);

    if (vector_size > 1) {

        auto **CHAM_descsubC11 = (CHAM_desc_t **) &pDescriptorC[1];
        auto **CHAM_descsubC12 = (CHAM_desc_t **) &pDescriptorC[2];
        auto **CHAM_descsubC22 = (CHAM_desc_t **) &pDescriptorC[3];

        *CHAM_descsubC11 = chameleon_desc_submatrix(*(CHAM_desc_t **) CHAM_descriptorC, 0, 0,
                                                    (*(CHAM_desc_t **) CHAM_descriptorC)->m / 2,
                                                    (*(CHAM_desc_t **) CHAM_descriptorC)->n / 2);
        *CHAM_descsubC12 = chameleon_desc_submatrix(*(CHAM_desc_t **) CHAM_descriptorC,
                                                    (*(CHAM_desc_t **) CHAM_descriptorC)->m / 2, 0,
                                                    (*(CHAM_desc_t **) CHAM_descriptorC)->m / 2,
                                                    (*(CHAM_desc_t **) CHAM_descriptorC)->n / 2);
        *CHAM_descsubC22 = chameleon_desc_submatrix(*(CHAM_desc_t **) CHAM_descriptorC,
                                                    (*(CHAM_desc_t **) CHAM_descriptorC)->m / 2,
                                                    (*(CHAM_desc_t **) CHAM_descriptorC)->n / 2,
                                                    (*(CHAM_desc_t **) CHAM_descriptorC)->m / 2,
                                                    (*(CHAM_desc_t **) CHAM_descriptorC)->n / 2);
    }
    ExageostatAllocateMatrixTile((&pDescriptorZ[0]), is_OOC, nullptr,
                                 (cham_flttype_t) float_point, dts, dts, dts * dts, N, 1, 0, 0, N, 1, p_grid,
                                 q_grid);
    ExageostatAllocateMatrixTile(pDescriptorZcpy, is_OOC, Zcpy, (cham_flttype_t) float_point, dts, dts,
                                 dts * dts, N, 1, 0, 0, N, 1, p_grid, q_grid);
    ExageostatAllocateMatrixTile(pDescriptorDeterminant, is_OOC, &dot_product_value, (cham_flttype_t) float_point,
                                 dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);

    for (int idx = 1; idx < pDescriptorZ.size(); idx++) {
        auto **CHAM_descriptorZ_ = &pDescriptorZ[idx];
        ExageostatAllocateMatrixTile(CHAM_descriptorZ_, is_OOC, nullptr, (cham_flttype_t) float_point, dts, dts,
                                     dts * dts, N / 2, 1, 0, 0, N / 2, 1, p_grid, q_grid);
    }

    for (auto &idx: pDescriptorProduct) {
        auto **CHAM_descriptorProduct = &idx;

       // EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descproduct, &data->dotp, ChamRealDouble, dts, dts, dts * dts, 1, 1, 0, 0, 1,
                                        //1, p_grid, q_grid);
        ExageostatAllocateMatrixTile(CHAM_descriptorProduct, is_OOC, &dot_product_value,
                                     (cham_flttype_t) float_point, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid,
                                     q_grid);
    }

    this->mpConfigurations->SetSequence(pSequence);
    this->mpConfigurations->SetRequest(request);

    //stop gsl error handler
    gsl_set_error_handler_off();
    free(Zcpy);
}

template<typename T>
void ChameleonImplementationDense<T>::ExaGeoStatInitContext(const int &apCoresNumber, const int &apGPUs) {

    if (!this->apContext) {
        CHAMELEON_user_tag_size(31, 26);
        CHAMELEON_Init(apCoresNumber, apGPUs)
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
void ChameleonImplementationDense<T>::CovarianceMatrixCodelet(void *apDescriptor, int &aTriangularPart,
                                                              dataunits::Locations *apLocation1,
                                                              dataunits::Locations *apLocation2,
                                                              dataunits::Locations *apLocation3,
                                                              double *aLocalTheta, int aDistanceMetric,
                                                              exageostat::kernels::Kernel *apKernel) {

    // Check for Initialise the Chameleon context.
    if (!this->apContext) {
        throw std::runtime_error(
                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStat<double/float>::ExaGeoStatInitializeHardware(configurations)'.");
    }

    RUNTIME_option_t options;
    RUNTIME_options_init(&options, (CHAM_context_t *) this->apContext,
                         (RUNTIME_sequence_t *) this->mpConfigurations->GetSequence(),
                         (RUNTIME_request_t *) this->mpConfigurations->GetRequest());

    int tempmm, tempnn;
    auto *CHAM_apDescriptor = (CHAM_desc_t *) apDescriptor;
    CHAM_desc_t A = *CHAM_apDescriptor;

    struct starpu_codelet *cl = &this->cl_dcmg;
    int m, n, m0 = 0, n0 = 0;

    cout << "Calling codlet" << endl;
    cout << "A.nt: " << A.nt << endl;
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
                               STARPU_VALUE, &apLocation1, sizeof(dataunits::Locations *),
                               STARPU_VALUE, &apLocation2, sizeof(dataunits::Locations *),
                               STARPU_VALUE, &apLocation3, sizeof(dataunits::Locations *),
                               STARPU_VALUE, &aLocalTheta, sizeof(double *),
                               STARPU_VALUE, &aDistanceMetric, sizeof(int),
                               STARPU_VALUE, &apKernel, sizeof(exageostat::kernels::Kernel *),
                               0);

            auto handle = (starpu_data_handle_t) RUNTIME_data_getaddr(CHAM_apDescriptor, m, n);
            this->apMatrix = (double *) starpu_variable_get_local_ptr(handle);
        }
    }
    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, (CHAM_context_t *) this->apContext);

    CHAMELEON_Sequence_Wait((RUNTIME_sequence_t *) this->mpConfigurations->GetSequence());

}

template<typename T>
void ChameleonImplementationDense<T>::GenerateObservationsVector(void *apDescriptor, Locations *apLocation1,
                                                                 Locations *apLocation2, Locations *apLocation3,
                                                                 vector<double> aLocalTheta, int aDistanceMetric,
                                                                 Kernel *apKernel) {

    // Check for Initialise the Chameleon context.
    if (!this->apContext) {
        throw std::runtime_error(
                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStat<double/float>::ExaGeoStatInitializeHardware(configurations)'.");
    }

    const int N = this->mpConfigurations->GetProblemSize();
    int seed = this->mpConfigurations->GetSeed();
    int iseed[4] = {seed, seed, seed, 1};

    //nomral random generation of e -- ei~N(0, 1) to generate Z
    auto *Nrand = (double *) malloc(N * sizeof(double));
    LAPACKE_dlarnv(3, iseed, N, Nrand);

    //Generate the co-variance matrix C
    auto *theta = (double *) malloc(aLocalTheta.size() * sizeof(double));
    for (int i = 0; i < aLocalTheta.size(); i++) {
        theta[i] = aLocalTheta[i];
    }

    VERBOSE("Initializing Covariance Matrix (Synthetic Dataset Generation Phase).....")
    int upper_lower = EXAGEOSTAT_LOWER;
    this->CovarianceMatrixCodelet(apDescriptor, upper_lower, apLocation1, apLocation2, apLocation3, theta,
                                  aDistanceMetric, apKernel);

    free(theta);
    VERBOSE("Done.\n")

    //Copy Nrand to Z
    VERBOSE("Generate Normal Random Distribution Vector Z (Synthetic Dataset Generation Phase) .....")
    auto **CHAM_descriptorZ = (CHAM_desc_t **) &this->mpConfigurations->GetDescriptorZ()[0];
    CopyDescriptorZ(*CHAM_descriptorZ, Nrand);
    VERBOSE("Done.\n")

    //Cholesky factorization for the Co-variance matrix C
    VERBOSE("Cholesky factorization of Sigma (Synthetic Dataset Generation Phase) .....")


    for (int i = 0; i < 10; i ++){
        auto matrix = ((CHAM_desc_t *) apDescriptor)->mat;
        cout << ((double *)matrix)[i] << endl;
    }


    int potential_failure = CHAMELEON_dpotrf_Tile(ChamLower, (CHAM_desc_t *) apDescriptor);
    FAILURE_LOGGER(potential_failure, "Factorization cannot be performed..\nThe matrix is not positive definite")
    VERBOSE("Done.\n")

    //Triangular matrix-matrix multiplication
    VERBOSE("Triangular matrix-matrix multiplication Z=L.e (Synthetic Dataset Generation Phase) .....")
    CHAMELEON_dtrmm_Tile(ChamLeft, ChamLower, ChamNoTrans, ChamNonUnit, 1, (CHAM_desc_t *) apDescriptor,
                         *CHAM_descriptorZ);
    VERBOSE("Done.\n")

    const int P = this->mpConfigurations->GetP();
    if (this->mpConfigurations->GetLogger()) {
        T *pMatrix;
        VERBOSE("Writing generated data to the disk (Synthetic Dataset Generation Phase) .....")
#ifdef CHAMELEON_USE_MPI
        pMatrix = (T*) malloc(N * sizeof(T));
        CHAMELEON_Tile_to_Lapack( *CHAM_descriptorZ, pMatrix, N);
        if ( CHAMELEON_My_Mpi_Rank() == 0 ){
            DiskWriter<T>::WriteVectorsToDisk(pMatrix, &N, &P, this->mpConfigurations->GetLoggerPath(), apLocation1);
        }
        free(pMatrix);
#else
        pMatrix = (T *) (*CHAM_descriptorZ)->mat;
        DiskWriter<T>::WriteVectorsToDisk(pMatrix, &N, &P, this->mpConfigurations->GetLoggerPath(), apLocation1);
        free(pMatrix);
#endif
        VERBOSE(" Done.\n")
    }

    CHAMELEON_dlaset_Tile(ChamUpperLower, 0, 0, (CHAM_desc_t *) apDescriptor);
    free(Nrand);
    VERBOSE("Done Z Vector Generation Phase. (Chameleon Synchronous)")
}

template<typename T>
void
ChameleonImplementationDense<T>::CopyDescriptorZ(void *apDescriptor, double *apDoubleVector) {

    // Check for Initialise the Chameleon context.
    if (!this->apContext) {
        throw std::runtime_error(
                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStat<double/float>::ExaGeoStatInitializeHardware(configurations)'.");
    }

    RUNTIME_option_t options;
    RUNTIME_options_init(&options, (CHAM_context_t *) this->apContext,
                         (RUNTIME_sequence_t *) this->mpConfigurations->GetSequence(),
                         (RUNTIME_request_t *) this->mpConfigurations->GetRequest());

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
void ChameleonImplementationDense<T>::DestoryDescriptors() {

    vector<void *> &pDescriptorC = this->mpConfigurations->GetDescriptorC();
    vector<void *> &pDescriptorZ = this->mpConfigurations->GetDescriptorZ();
    auto pChameleonDescriptorZcpy = (CHAM_desc_t **) &this->mpConfigurations->GetDescriptorZcpy();
    vector<void *> &pDescriptorProduct = this->mpConfigurations->GetDescriptorProduct();
    auto pChameleonDescriptorDeterminant = (CHAM_desc_t **) &this->mpConfigurations->GetDescriptorDeterminant();

    if (!pDescriptorC.empty() && pDescriptorC[0]) {
        CHAMELEON_Desc_Destroy((CHAM_desc_t **) &pDescriptorC[0]);
    }
    if (!pDescriptorZ.empty() && pDescriptorZ[0]) {
        CHAMELEON_Desc_Destroy((CHAM_desc_t **) &pDescriptorZ[0]);
    }
    if (!pDescriptorProduct.empty() && pDescriptorProduct[0]) {
        CHAMELEON_Desc_Destroy((CHAM_desc_t **) &pDescriptorProduct[0]);
    }
    if (*pChameleonDescriptorZcpy) {
        CHAMELEON_Desc_Destroy(pChameleonDescriptorZcpy);
    }
    if (*pChameleonDescriptorDeterminant) {
        CHAMELEON_Desc_Destroy(pChameleonDescriptorDeterminant);
    }

    if ((RUNTIME_sequence_t *) this->mpConfigurations->GetSequence()) {
        CHAMELEON_Sequence_Destroy((RUNTIME_sequence_t *) this->mpConfigurations->GetSequence());
    }
}

template<typename T>
void ChameleonImplementationDense<T>::ExageostatAllocateMatrixTile(void **apDescriptor, bool ais_OOC, T *apMemSpace,
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
T ChameleonImplementationDense<T>::ExageostatMleTile(vector<double> &apTheta, Locations * apDataLocations, DataModelingConfigurations *apDataModelingConfiguration){

    //Initialization
    auto median_locations = apDataLocations->CalculateMedianLocations();
    T loglik = 0.0, logdet = 0.0, time_facto = 0.0, time_solve = 0.0, logdet_calculate = 0.0, matrix_gen_time = 0.0;
    T dzcpy_time = 0.0;

    int n, nhrs, success, i, num_params;
    T flops = 0.0;
    T*univariate_theta;
    T*univariate2_theta;
    T*univariate3_theta;
    T nu12;
    T rho;
    T sigma_square12;
    *apDataModelingConfiguration->GetDeterminantValue() = 0;

    vector<void*> &descriptor_c =  this->mpConfigurations->GetDescriptorC();
    auto **EXA_descsubC11 = &descriptor_c[1];
    auto **EXA_descsubC12 =  &descriptor_c[2];
    auto **EXA_descsubC22 = &descriptor_c[3];
    vector<void*> &descriptor_z =  this->mpConfigurations->GetDescriptorZ();
    auto **EXA_descZ1 =  &descriptor_c[1];
    auto **EXA_descZ2 =   &descriptor_c[2];
    auto pDescriptorZcpy = &this->mpConfigurations->GetDescriptorZcpy();

    auto EXA_descdet =  &this->mpConfigurations->GetDescriptorDeterminant();
    vector<void*> &descriptor_product =  this->mpConfigurations->GetDescriptorProduct();
    void *EXA_descproduct1  =   this->mpConfigurations->GetDescriptorProduct()[1];
    void *EXA_descproduct2 =  this->mpConfigurations->GetDescriptorProduct()[2];
    auto msequence =  (RUNTIME_sequence_t *)this->mpConfigurations->GetSequence();
    auto mrequest = (RUNTIME_request_t *) this->mpConfigurations->GetRequest();
    void *EXA_descriptorC = this->mpConfigurations->GetDescriptorC()[0];
    void *EXA_descriptorZ = this->mpConfigurations->GetDescriptorZ()[0];
    void *EXA_descZcpy =  this->mpConfigurations->GetDescriptorZcpy();
    void *EXA_descriptor_product = this->mpConfigurations->GetDescriptorProduct()[0];

    num_params = apDataModelingConfiguration->GetParametersNumber();
    //TODO: make custom getters for both n and nhrs depending on linear algebra solver
    cout << "after get param" <<endl;

    n = ((CHAM_desc_t *) EXA_descriptorC)->m;
    nhrs = (*(CHAM_desc_t **) EXA_descriptorZ)->n;

    auto test = (CHAM_desc_t *) EXA_descriptorZ;

    START_TIMING(dzcpy_time);
    int iter_count = apDataModelingConfiguration->GetIterationsValue();
    if (iter_count == 0) {
        //Save a copy of descZ into descZcpy for restoring each iteration (Only for the first iteration)
        ExageostatLacpyTile(UpperLower::EXAGEOSTAT_UPPER_LOWER, EXA_descriptorZ, EXA_descZcpy);
    }
    string recovery_file = apDataModelingConfiguration->GetRecoveryFile();

    if ( recovery_file.empty() || !(this->recover((char *)(recovery_file.c_str()), iter_count, (T*) apTheta.data(), &loglik, num_params))) {
        START_TIMING(dzcpy_time);
        if (iter_count == 0) {
            //Save a copy of descZ into descZcpy for restoring each iteration (Only for the first iteration)
            ExageostatLacpyTile(EXAGEOSTAT_UPPER_LOWER, EXA_descriptorZ, EXA_descZcpy);
        } else {
            VERBOSE("Re-store the original Z vector...");
            ExageostatLacpyTile(EXAGEOSTAT_UPPER_LOWER, EXA_descZcpy, EXA_descriptorZ);
            VERBOSE(" Done.\n")
        }
        STOP_TIMING(dzcpy_time)

        //Generate new co-variance matrix C based on new theta
        VERBOSE("Generate New Covariance Matrix...")
        START_TIMING(matrix_gen_time)
        auto kernel_fun = apDataModelingConfiguration->GetKernel();

        if (kernel_fun == "bivariate_matern_parsimonious2" ||
            kernel_fun == "bivariate_matern_parsimonious2_profile") {

            int upper_lower = EXAGEOSTAT_UPPER_LOWER;
            exageostat::kernels::Kernel * kernel = exageostat::plugins::PluginRegistry<exageostat::kernels::Kernel>::Create(
                    kernel_fun);

            univariate_theta = (T *) malloc(3 * sizeof(T));
            univariate2_theta = (T *) malloc(3 * sizeof(T));
            univariate3_theta = (T *) malloc(3 * sizeof(T));
            univariate_theta[0] = apTheta[0];
            univariate_theta[1] = apTheta[2];
            univariate_theta[2] = apTheta[3];

            std::string distance_metric = apDataModelingConfiguration->GetDistanceMetric();
            std::string kernel_name = "univariate_matern_stationary";

            this->CovarianceMatrixCodelet(EXA_descsubC11, upper_lower, apDataLocations, apDataLocations,
                                          median_locations, (double*)univariate_theta,
                                          0, kernel);

            nu12 = 0.5 * (apTheta[3] + apTheta[4]);

            rho = apTheta[5] * sqrt((tgamma(apTheta[3] + 1) * tgamma(apTheta[4] + 1)) /
                                    (tgamma(apTheta[3]) * tgamma(apTheta[4]))) *
                  tgamma(nu12) / tgamma(nu12 + 1);
            sigma_square12 = rho * sqrt(apTheta[0] * apTheta[1]);

            univariate2_theta[0] = sigma_square12;
            univariate2_theta[1] = apTheta[2];
            univariate2_theta[2] = nu12;

            this->CovarianceMatrixCodelet(EXA_descsubC12, upper_lower, median_locations, apDataLocations,
                                          median_locations, (double*)univariate2_theta,
                                          0, kernel);

            STOP_TIMING(matrix_gen_time);
            VERBOSE(" Done.\n");

            univariate3_theta[0] = apTheta[1];
            univariate3_theta[1] = apTheta[2];
            univariate3_theta[2] = apTheta[4];

            this->CovarianceMatrixCodelet(EXA_descsubC22, upper_lower, median_locations, apDataLocations,
                                          median_locations, (double*)univariate2_theta,
                                          0, kernel);
        } else {
            std::string distance_metric = apDataModelingConfiguration->GetDistanceMetric();
            std::string kernel_name = apDataModelingConfiguration->GetKernel();
            cout << "K: " << kernel_name << endl;


            int upper_lower = EXAGEOSTAT_LOWER;

            cout << "kernel_fun: " << kernel_fun << endl;
            exageostat::kernels::Kernel * kernel = exageostat::plugins::PluginRegistry<exageostat::kernels::Kernel>::Create(
                    kernel_name);

            cout << " Inside here" << endl;
            cout << kernel->GetPValue() << endl;
            cout << "mle theta [0] = "<< apTheta[0] << endl;
            this->CovarianceMatrixCodelet(EXA_descriptorC, upper_lower, apDataLocations, apDataLocations,
                                          apDataLocations, apTheta.data(),
                                          0, kernel);
        }
        ExageostatSequenceWait(msequence);
        STOP_TIMING(matrix_gen_time);
        VERBOSE(" Done.\n");

        VERBOSE("Cholesky factorization of Sigma...");
        START_TIMING(time_facto);

        for (int i = 0; i < 10; i ++){
            auto matrix = ((CHAM_desc_t *) EXA_descriptorC)->mat;
            cout << ((double *)matrix)[i] << endl;
        }

        success = CHAMELEON_dpotrf_Tile(ChamLower, (CHAM_desc_t *) EXA_descriptorC);
        STOP_TIMING(time_facto);
        FAILURE_LOGGER(success, "Factorization cannot be performed..\n The matrix is not positive definite\n\n");
        flops = flops + FLOPS_DPOTRF(n);
        VERBOSE(" Done.\n");

        //Calculate log(|C|) --> log(square(|L|))
        VERBOSE("Calculating the log determinant ...");
        START_TIMING(logdet_calculate);
        ExageostatMleMdetTileAsync(EXA_descriptorC, msequence, &mrequest[0], *EXA_descdet);
        ExageostatSequenceWait(msequence);

        // printf("det: %f\n", data->det);
        logdet = 2 * (*apDataModelingConfiguration->GetDeterminantValue());
        STOP_TIMING(logdet_calculate);
        VERBOSE(" Done.\n");

        //Solving Linear System (L*X=Z)--->inv(L)*Z
        VERBOSE("Solving the linear system ...\n");
        START_TIMING(time_solve);
        ExageostatDtrsmTile(EXAGEOSTAT_LEFT, EXAGEOSTAT_LOWER, EXAGEOSTAT_NO_TRANS, EXAGEOSTAT_NON_UNIT, 1, EXA_descriptorC, EXA_descriptorZ);
        STOP_TIMING(time_solve);
        flops = flops + FLOPS_DTRSM(ChamLeft, n, nhrs);
        VERBOSE(" Done.\n")

        //Calculate MLE likelihood
        VERBOSE("Calculating the MLE likelihood function ...");

        ExageostatDgemmTile(EXAGEOSTAT_TRANS, EXAGEOSTAT_NO_TRANS, 1, EXA_descriptorZ, EXA_descriptorZ, 0, EXA_descriptor_product);

        if (apDataModelingConfiguration->GetKernel() == "bivariate_matern_parsimonious2_profile") {

            loglik = -(n / 2) + (n / 2) * log(n) - (n / 2) * log(*apDataModelingConfiguration->GetDotProduct()[0]) - 0.5 * logdet -
                     (T) (n / 2.0) * log(2.0 * PI);
            ExageostatDgemmTile(EXAGEOSTAT_TRANS, EXAGEOSTAT_NO_TRANS, 1, EXA_descZ1, EXA_descZ1, 0, EXA_descproduct1);
            ExageostatDgemmTile(EXAGEOSTAT_TRANS, EXAGEOSTAT_NO_TRANS, 1, EXA_descZ2, EXA_descZ2, 0, EXA_descproduct2);

            ///where do dotp1 and 2 get initialized???
            *apDataModelingConfiguration->GetVariance()[1] = (1.0 / (n / 2)) * (*apDataModelingConfiguration->GetDotProduct()[1]);
            *apDataModelingConfiguration->GetVariance()[2]  = (1.0 / (n / 2)) * (*apDataModelingConfiguration->GetDotProduct()[2]);

        }else {
            loglik = -0.5 * (*apDataModelingConfiguration->GetDotProduct()[1]) - 0.5 * logdet - (double) (n / 2.0) * log(2.0 * PI);
            *(apDataModelingConfiguration->GetVariance()[0]) = apTheta[0];
        }
        VERBOSE(" Done.\n");
    }

    fprintf(stderr, " %3d- Model Parameters (", iter_count + 1);
    if (apDataModelingConfiguration->GetLog())
        ///should we be leaving the file opening and closing to the user?
        fprintf(apDataModelingConfiguration->GetFileLog(), " %3d- Model Parameters (",iter_count + 1);

    if ((apDataModelingConfiguration->GetKernel() == "bivariate_matern_parsimonious_profile") ||
        (apDataModelingConfiguration->GetKernel() == "bivariate_matern_parsimonious2_profile")) {
        fprintf(stderr, "%.8f, %.8f,", *apDataModelingConfiguration->GetVariance()[1], *apDataModelingConfiguration->GetVariance()[2]);
        if (apDataModelingConfiguration->GetLog())
            fprintf(apDataModelingConfiguration->GetFileLog(), "%.8f, %.8f,", *apDataModelingConfiguration->GetVariance()[1], *apDataModelingConfiguration->GetVariance()[2]);
        i = 2;
    } else {
        i = 0;
    }
    for (; i < num_params; i++) {
        fprintf(stderr, "%.8f", apTheta[i]);
        if (i < num_params - 1) {
            fprintf(stderr, ",");
        }

        if (apDataModelingConfiguration->GetLog()) {
            fprintf(apDataModelingConfiguration->GetFileLog(), "%.8f, ", apTheta[i]);
        }
    }
    fprintf(stderr, ")----> LogLi: %.18f\n", loglik);
    if (apDataModelingConfiguration->GetLog()) {
        fprintf(apDataModelingConfiguration->GetFileLog(), ")----> LogLi: %.18f\n", loglik);
    }


    fprintf(stderr, " ---- Facto Time: %6.2f seconds\n", time_facto);
    fprintf(stderr, " ---- Matrix Generation Time: %6.2f\n", matrix_gen_time);
    fprintf(stderr, " ---- Total Time: %6.2f\n", matrix_gen_time + time_facto + logdet_calculate + time_solve);

    apDataModelingConfiguration->SetIterationsValue(apDataModelingConfiguration->GetIterationsValue() + 1);
    // for experiments
    if(DataModelingConfigurations::GetRunMode() == VERBOSE_MODE) {
        apDataModelingConfiguration->SetAvgExecutedTimePerIteration( apDataModelingConfiguration->GetAvgExecutedTimePerIteration() +/*matrix_gen_time*/+time_facto + logdet_calculate + time_solve);
        apDataModelingConfiguration->SetAvgFlopsPerIteration(apDataModelingConfiguration->GetAvgFlopsPerIteration() + flops / 1e9 / (time_facto + time_solve));
    }
     apDataModelingConfiguration->SetFinalLogLik(loglik);

    //output
////    results.final_loglik = loglik;
    return loglik;
}

template<typename T>
int ChameleonImplementationDense<T>::ExageostatLacpyTile(exageostat::common::UpperLower aUpperLower, void *apA, void *apB){
    cout << "LPCY  descriptor m: " << ((CHAM_desc_t *) apA)->m << " n: " << ((CHAM_desc_t *) apA)->n << endl;
    cout << "LPCY B descriptor m: " << ((CHAM_desc_t *) apB)->m << " n: " << ((CHAM_desc_t *) apB)->n << endl;

    return CHAMELEON_dlacpy_Tile((cham_uplo_t) aUpperLower, (CHAM_desc_t *) apA, (CHAM_desc_t *) apB);
}

template<typename T>
int ChameleonImplementationDense<T>::ExageostatLap2Desc(exageostat::common::UpperLower aUpperLower, void *apAf77, int aLda, void * apA){
    return CHAMELEON_Lap2Desc((cham_uplo_t) aUpperLower, (CHAM_desc_t *) apAf77, aLda, (CHAM_desc_t *) apA);
}


template<typename T>
int ChameleonImplementationDense<T>::ExageostatSequenceWait(void * apSequence) {
    return CHAMELEON_Sequence_Wait((RUNTIME_sequence_t *)apSequence);
}

template<typename T>
int ChameleonImplementationDense<T>::ExageostatDpotrfTile(exageostat::common::UpperLower aUpperLower, void *apA){
    return CHAMELEON_dpotrf_Tile((cham_uplo_t) aUpperLower, (CHAM_desc_t *) apA);
}

template<typename T>
int ChameleonImplementationDense<T>::ExageostatDtrsmTile(common::Side aSide, common::UpperLower aUpperLower, common::Trans aTrans, common::Diag aDiag, T aAlpha, void *apA, void *apB) {
    return CHAMELEON_dtrsm_Tile((cham_side_t) aSide, (cham_uplo_t) aUpperLower, (cham_trans_t) aTrans, (cham_diag_t) aDiag, aAlpha, (CHAM_desc_t *) apA, (CHAM_desc_t *) apB);
}

template<typename T>
int ChameleonImplementationDense<T>::ExageostatDgemmTile(common::Trans aTransA, common::Trans aTransB, T aAlpha,
                                                         void *apA, void *apB, T aBeta, void *apC) {
    return CHAMELEON_dgemm_Tile((cham_trans_t) aTransA, (cham_trans_t) aTransB, aAlpha,(CHAM_desc_t *) apA, (CHAM_desc_t *) apB, aBeta, (CHAM_desc_t *)apC);
}

template<typename T>
int ChameleonImplementationDense<T>::ExageostatMleMdetTileAsync(void *apDescA, void * apSequence, void *apRequest, void *apDescDet){
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    chamctxt = chameleon_context_self();

    if (((RUNTIME_sequence_t *)apSequence)->status != CHAMELEON_SUCCESS) {
        return -2;
    }
    RUNTIME_options_init(&options, chamctxt, (RUNTIME_sequence_t *)apSequence, (RUNTIME_request_t *)apRequest);
    int m, m0, n0;
    int tempmm;
    CHAM_desc_t A = *((CHAM_desc_t *) apDescA);
    struct starpu_codelet *cl = &this->cl_dmdet;

    for (m = 0; m < A.mt; m++) {
        tempmm = m == A.mt - 1 ? A.m - m * A.mb : A.mb;
        starpu_insert_task(cl,
                           STARPU_VALUE, &tempmm, sizeof(int),
                           STARPU_VALUE, &tempmm, sizeof(int),
                           STARPU_R, (starpu_data_handle_t)RUNTIME_data_getaddr( (CHAM_desc_t *)apDescA, m, 0),
                           STARPU_VALUE, &m0, sizeof(int),
                           STARPU_VALUE, &n0, sizeof(int),
                           STARPU_RW, (starpu_data_handle_t)RUNTIME_data_getaddr( (CHAM_desc_t *)apDescDet, 0, 0),
                           0);
    }

    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, chamctxt);
    return CHAMELEON_SUCCESS;
}


namespace exageostat::linearAlgebra::dense {
    template<typename T> void *ChameleonImplementationDense<T>::apContext = nullptr;
}



template<typename T>
int ChameleonImplementationDense<T>::ExageostaStrideVecTileAsync(void *apDescA, void *apDescB, void *apDescC,
                                                                 void * apSequence, void *apRequest){
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    chamctxt = chameleon_context_self();
    if (((RUNTIME_sequence_t *)apSequence)->status != CHAMELEON_SUCCESS)
        return -2;
    RUNTIME_options_init(&options, chamctxt, (RUNTIME_sequence_t *)apSequence, (RUNTIME_request_t *)apRequest);

    int m, m0;
    int tempmm;
    CHAM_desc_t A = *((CHAM_desc_t *)apDescA);
    CHAM_desc_t B = *((CHAM_desc_t *)apDescB);
    CHAM_desc_t C = *((CHAM_desc_t *)apDescC);
    struct starpu_codelet *cl = &this->cl_stride_vec;
    for (m = 0; m < A.mt; m++) {
        tempmm = m == A.mt - 1 ? A.m - m * A.mb : A.mb;
        m0 = m * A.mb;
        starpu_insert_task(cl,
                           STARPU_VALUE, &tempmm, sizeof(int),
                           STARPU_VALUE, &m0, sizeof(int),
                           STARPU_VALUE, &m, sizeof(int),
                           STARPU_R, (starpu_data_handle_t)RUNTIME_data_getaddr((CHAM_desc_t *)apDescA, m, 0 ),
                           STARPU_W, (starpu_data_handle_t)RUNTIME_data_getaddr((CHAM_desc_t *)apDescB, (int) floor(m / 2.0), 0 ),
                           STARPU_W, (starpu_data_handle_t)RUNTIME_data_getaddr((CHAM_desc_t *)apDescC, (int) floor(m / 2.0), 0 ),
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
                STARPU_NAME, "stride_vec",
#endif
                           0);

    }
    RUNTIME_options_ws_free(&options);
    return CHAMELEON_SUCCESS;
}


template<typename T>
int ChameleonImplementationDense<T>::ExageostatMleDcmgTileAsync(common::UpperLower aUpperLower, void *apDescA, dataunits::Locations *apL1,
                                                                dataunits::Locations *apL2, dataunits::Locations *apLm, T* apTheta,
                                                                std::string &aDm, std::string &aKernelFun,
                                                                void *apSequence, void *apRequest){
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    chamctxt = chameleon_context_self();

    if (((RUNTIME_sequence_t *)apSequence)->status != CHAMELEON_SUCCESS) {
        cout << "dcmg tile async returened -2\n";
        return -2; }
    RUNTIME_options_init(&options, chamctxt, ((RUNTIME_sequence_t *)apSequence), ((RUNTIME_request_t *)apRequest));

    int m, n, m0, n0;
    int distance_metric = aDm == "gc" ? 1 : 0;
    int kernel;
    int tempmm, tempnn;
    CHAM_desc_t A = *((CHAM_desc_t *)apDescA);
    struct starpu_codelet *cl = &this->cl_dcmg;
    int size = A.n;

    for (n = 0; n < A.nt; n++) {
        tempnn = n == A.nt - 1 ? A.n - n * A.nb : A.nb;
        if (aUpperLower == ChamUpperLower)
            m = 0;
        else
            m = A.m == A.n ? n : 0;
        for (; m < A.mt; m++) {

            tempmm = m == A.mt - 1 ? A.m - m * A.mb : A.mb;
            m0 = m * A.mb;
            n0 = n * A.nb;
            starpu_insert_task(cl,
                               STARPU_VALUE, &tempmm, sizeof(int),
                               STARPU_VALUE, &tempnn, sizeof(int),
                               STARPU_VALUE, &m0, sizeof(int),
                               STARPU_VALUE, &n0, sizeof(int),
                               STARPU_W, (starpu_data_handle_t)RUNTIME_data_getaddr((CHAM_desc_t *)apDescA, m, n ),
                               STARPU_VALUE, &apL1, sizeof(Locations * ),
                               STARPU_VALUE, &apL2, sizeof(Locations * ),
                               STARPU_VALUE, &apLm, sizeof(Locations * ),
                               STARPU_VALUE, &apTheta, sizeof(double* ),
                               STARPU_VALUE, &distance_metric, sizeof(int),
                               STARPU_VALUE, &aKernelFun, sizeof(int),
                               0);
        }

    }
}
