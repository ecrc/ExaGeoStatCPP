
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

using namespace std;

using namespace exageostat::linearAlgebra::dense;
using namespace exageostat::common;
using namespace exageostat::dataunits;
using namespace exageostat::kernels;
using namespace exageostat::helpers;
using namespace exageostat::configurations;

// Define a method to set up the Chameleon descriptors
template<typename T>
void ChameleonImplementationDense<T>::InitiateDescriptors() {

    // Check for Initialise the Chameleon context.
    if (!this->apContext) {
        throw std::runtime_error(
                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStat<double/float>::ExaGeoStatInitializeHardware(configurations)'.");
    }

    // Get the common Dictionary of the Configurations
    configurations::Configurations *configurations = configurations::Configurations::GetConfigurations();
    
    // Declare variables for Chameleon descriptors
    vector<void *> &pDescriptorC = configurations->GetDescriptorC();
    vector<void *> &pDescriptorZ = configurations->GetDescriptorZ();
    auto pDescriptorZcpy = &configurations->GetDescriptorZcpy();
    vector<void *> &pDescriptorProduct = configurations->GetDescriptorProduct();
    auto pDescriptorDeterminant = &configurations->GetDescriptorDeterminant();

    // Get the problem size and other configuration parameters
    int vector_size;
    int N = configurations->GetProblemSize();
    int dts = configurations->GetDenseTileSize();
    int p_grid = configurations->GetPGrid();
    int q_grid = configurations->GetQGrid();
    bool is_OOC = configurations->GetIsOOC();

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
    ExaGeoStatAllocateMatrixTile(CHAM_descriptorC, is_OOC, nullptr, (cham_flttype_t) float_point, dts, dts,
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
    ExaGeoStatAllocateMatrixTile((&pDescriptorZ[0]), is_OOC, nullptr,
                                 (cham_flttype_t) float_point, dts, dts, dts * dts, N, 1, 0, 0, N, 1, p_grid,
                                 q_grid);
    ExaGeoStatAllocateMatrixTile(pDescriptorZcpy, is_OOC, Zcpy, (cham_flttype_t) float_point, dts, dts,
                                 dts * dts, N, 1, 0, 0, N, 1, p_grid, q_grid);
    ExaGeoStatAllocateMatrixTile(pDescriptorDeterminant, is_OOC, &dot_product_value, (cham_flttype_t) float_point,
                                 dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid, q_grid);

    for (int idx = 1; idx < pDescriptorZ.size(); idx++) {
        auto **CHAM_descriptorZ_ = &pDescriptorZ[idx];
        ExaGeoStatAllocateMatrixTile(CHAM_descriptorZ_, is_OOC, nullptr, (cham_flttype_t) float_point, dts, dts,
                                     dts * dts, N / 2, 1, 0, 0, N / 2, 1, p_grid, q_grid);
    }

    for (auto &idx: pDescriptorProduct) {
        auto **CHAM_descriptorProduct = &idx;

       // EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descproduct, &data->dotp, ChamRealDouble, dts, dts, dts * dts, 1, 1, 0, 0, 1,
                                        //1, p_grid, q_grid);
        ExaGeoStatAllocateMatrixTile(CHAM_descriptorProduct, is_OOC, &dot_product_value,
                                     (cham_flttype_t) float_point, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, p_grid,
                                     q_grid);
    }

    configurations->SetSequence(pSequence);
    configurations->SetRequest(request);

    //stop gsl error handler
    gsl_set_error_handler_off();
    free(Zcpy);
}

template<typename T>
void ChameleonImplementationDense<T>::ExaGeoStatInitContext() {

    if (!this->apContext) {
        CHAMELEON_user_tag_size(31, 26);
        CHAMELEON_Init(Configurations::GetConfigurations()->GetCoresNumber(), Configurations::GetConfigurations()->GetGPUsNumbers())
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

    configurations::Configurations *configurations = configurations::Configurations::GetConfigurations();

    RUNTIME_option_t options;
    RUNTIME_options_init(&options, (CHAM_context_t *) this->apContext,
                         (RUNTIME_sequence_t *) configurations->GetSequence(),
                         (RUNTIME_request_t *) configurations->GetRequest());

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

    CHAMELEON_Sequence_Wait((RUNTIME_sequence_t *) configurations->GetSequence());

}

template<typename T>
void ChameleonImplementationDense<T>::GenerateObservationsVector(void *apDescriptor, Locations *apLocation1,
                                                                 Locations *apLocation2, Locations *apLocation3,
                                                                 vector<double> aLocalTheta, int aDistanceMetric,
                                                                 Kernel *apKernel) {

    configurations::Configurations *configurations = configurations::Configurations::GetConfigurations();

    // Check for Initialise the Chameleon context.
    if (!this->apContext) {
        throw std::runtime_error(
                "ExaGeoStat hardware is not initialized, please use 'ExaGeoStat<double/float>::ExaGeoStatInitializeHardware(configurations)'.");
    }

    const int N = configurations->GetProblemSize();
    int seed = configurations->GetSeed();
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
    auto **CHAM_descriptorZ = (CHAM_desc_t **) &configurations->GetDescriptorZ()[0];
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

    const int P = configurations->GetP();
    if (configurations->GetLogger()) {
        T *pMatrix;
        VERBOSE("Writing generated data to the disk (Synthetic Dataset Generation Phase) .....")
#ifdef CHAMELEON_USE_MPI
        pMatrix = (T*) malloc(N * sizeof(T));
        CHAMELEON_Tile_to_Lapack( *CHAM_descriptorZ, pMatrix, N);
        if ( CHAMELEON_My_Mpi_Rank() == 0 ){
            DiskWriter<T>::WriteVectorsToDisk(pMatrix, &N, &P, configurations->GetLoggerPath(), apLocation1);
        }
        free(pMatrix);
#else
        pMatrix = (T *) (*CHAM_descriptorZ)->mat;
        DiskWriter<T>::WriteVectorsToDisk(pMatrix, &N, &P, configurations->GetLoggerPath(), apLocation1);
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

    configurations::Configurations *configurations = configurations::Configurations::GetConfigurations();

    RUNTIME_option_t options;
    RUNTIME_options_init(&options, (CHAM_context_t *) this->apContext,
                         (RUNTIME_sequence_t *) configurations->GetSequence(),
                         (RUNTIME_request_t *) configurations->GetRequest());

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
void ChameleonImplementationDense<T>::DestroyDescriptors() {

    configurations::Configurations *configurations = configurations::Configurations::GetConfigurations();

    vector<void *> &pDescriptorC = configurations->GetDescriptorC();
    vector<void *> &pDescriptorZ = configurations->GetDescriptorZ();
    auto pChameleonDescriptorZcpy = (CHAM_desc_t **) &configurations->GetDescriptorZcpy();
    vector<void *> &pDescriptorProduct = configurations->GetDescriptorProduct();
    auto pChameleonDescriptorDeterminant = (CHAM_desc_t **) &configurations->GetDescriptorDeterminant();

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

    if ((RUNTIME_sequence_t *) configurations->GetSequence()) {
        CHAMELEON_Sequence_Destroy((RUNTIME_sequence_t *) configurations->GetSequence());
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
T ChameleonImplementationDense<T>::ExaGeoStatMleTile(Locations * apDataLocations){

    //Initialization
    configurations::Configurations *configurations = configurations::Configurations::GetConfigurations();

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
    //// TODO: Fix this!
//    *configurations->GetDeterminantValue() = 0;

    vector<void*> &descriptor_c =  configurations->GetDescriptorC();
    auto **EXA_descsubC11 = &descriptor_c[1];
    auto **EXA_descsubC12 =  &descriptor_c[2];
    auto **EXA_descsubC22 = &descriptor_c[3];
    vector<void*> &descriptor_z =  configurations->GetDescriptorZ();
    auto **EXA_descZ1 =  &descriptor_c[1];
    auto **EXA_descZ2 =   &descriptor_c[2];
    auto pDescriptorZcpy = &configurations->GetDescriptorZcpy();

    auto EXA_descdet =  &configurations->GetDescriptorDeterminant();
    vector<void*> &descriptor_product =  configurations->GetDescriptorProduct();
    void *EXA_descproduct1  =   configurations->GetDescriptorProduct()[1];
    void *EXA_descproduct2 =  configurations->GetDescriptorProduct()[2];
    auto msequence =  (RUNTIME_sequence_t *)configurations->GetSequence();
    auto mrequest = (RUNTIME_request_t *) configurations->GetRequest();
    void *EXA_descriptorC = configurations->GetDescriptorC()[0];
    void *EXA_descriptorZ = configurations->GetDescriptorZ()[0];
    void *EXA_descZcpy =  configurations->GetDescriptorZcpy();
    void *EXA_descriptor_product = configurations->GetDescriptorProduct()[0];
    auto theta = configurations->GetInitialTheta();

    num_params = configurations->GetParametersNumber();
    //TODO: make custom getters for both n and nhrs depending on linear algebra solver
    cout << "after get param" <<endl;

    n = ((CHAM_desc_t *) EXA_descriptorC)->m;
    nhrs = (*(CHAM_desc_t **) EXA_descriptorZ)->n;

    auto test = (CHAM_desc_t *) EXA_descriptorZ;

    START_TIMING(dzcpy_time);
    int iter_count = configurations->GetIterationsValue();
    if (iter_count == 0) {
        //Save a copy of descZ into descZcpy for restoring each iteration (Only for the first iteration)
        ExaGeoStatLapackCopyTile(UpperLower::EXAGEOSTAT_UPPER_LOWER, EXA_descriptorZ, EXA_descZcpy);
    }
    string recovery_file = configurations->GetRecoveryFile();

    if ( recovery_file.empty() || !(this->recover((char *)(recovery_file.c_str()), iter_count, (T*) theta.data(), &loglik, num_params))) {
        START_TIMING(dzcpy_time);
        if (iter_count == 0) {
            //Save a copy of descZ into descZcpy for restoring each iteration (Only for the first iteration)
            ExaGeoStatLapackCopyTile(EXAGEOSTAT_UPPER_LOWER, EXA_descriptorZ, EXA_descZcpy);
        } else {
            VERBOSE("Re-store the original Z vector...");
            ExaGeoStatLapackCopyTile(EXAGEOSTAT_UPPER_LOWER, EXA_descZcpy, EXA_descriptorZ);
            VERBOSE(" Done.\n")
        }
        STOP_TIMING(dzcpy_time)

        //Generate new co-variance matrix C based on new theta
        VERBOSE("Generate New Covariance Matrix...")
        START_TIMING(matrix_gen_time)
        auto kernel_fun = configurations->GetKernelName();

        if (kernel_fun == "bivariate_matern_parsimonious2" ||
            kernel_fun == "bivariate_matern_parsimonious2_profile") {

            int upper_lower = EXAGEOSTAT_UPPER_LOWER;
            exageostat::kernels::Kernel * kernel = exageostat::plugins::PluginRegistry<exageostat::kernels::Kernel>::Create(
                    kernel_fun);

            univariate_theta = (T *) malloc(3 * sizeof(T));
            univariate2_theta = (T *) malloc(3 * sizeof(T));
            univariate3_theta = (T *) malloc(3 * sizeof(T));
            univariate_theta[0] = theta[0];
            univariate_theta[1] = theta[2];
            univariate_theta[2] = theta[3];

            //// TODO:: distance matrix
//            std::string distance_metric = configurations.GetDistanceMetric();
            std::string kernel_name = "univariate_matern_stationary";

            this->CovarianceMatrixCodelet(EXA_descsubC11, upper_lower, apDataLocations, apDataLocations,
                                          median_locations, (double*)univariate_theta,
                                          0, kernel);

            nu12 = 0.5 * (theta[3] + theta[4]);

            rho = theta[5] * sqrt((tgamma(theta[3] + 1) * tgamma(theta[4] + 1)) /
                                    (tgamma(theta[3]) * tgamma(theta[4]))) *
                  tgamma(nu12) / tgamma(nu12 + 1);
            sigma_square12 = rho * sqrt(theta[0] * theta[1]);

            univariate2_theta[0] = sigma_square12;
            univariate2_theta[1] = theta[2];
            univariate2_theta[2] = nu12;

            this->CovarianceMatrixCodelet(EXA_descsubC12, upper_lower, median_locations, apDataLocations,
                                          median_locations, (double*)univariate2_theta,
                                          0, kernel);

            STOP_TIMING(matrix_gen_time);
            VERBOSE(" Done.\n");

            univariate3_theta[0] = theta[1];
            univariate3_theta[1] = theta[2];
            univariate3_theta[2] = theta[4];

            this->CovarianceMatrixCodelet(EXA_descsubC22, upper_lower, median_locations, apDataLocations,
                                          median_locations, (double*)univariate2_theta,
                                          0, kernel);
        } else {
            std::string kernel_name = configurations->GetKernelName();
            cout << "K: " << kernel_name << endl;


            int upper_lower = EXAGEOSTAT_LOWER;

            cout << "kernel_fun: " << kernel_fun << endl;
            exageostat::kernels::Kernel * kernel = exageostat::plugins::PluginRegistry<exageostat::kernels::Kernel>::Create(
                    kernel_name);

            cout << " Inside here" << endl;
            cout << kernel->GetPValue() << endl;
            cout << "mle theta [0] = "<< theta[0] << endl;
            this->CovarianceMatrixCodelet(EXA_descriptorC, upper_lower, apDataLocations, apDataLocations,
                                          apDataLocations, theta.data(),
                                          0, kernel);
        }
        ExaGeoStatSequenceWait(msequence);
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
        ExaGeoStatMeasureDetTileAsync(EXA_descriptorC, msequence, &mrequest[0], *EXA_descdet);
        ExaGeoStatSequenceWait(msequence);

        // printf("det: %f\n", data->det);
        logdet = 2 * configurations->GetDeterminantValue();
        STOP_TIMING(logdet_calculate);
        VERBOSE(" Done.\n");

        //Solving Linear System (L*X=Z)--->inv(L)*Z
        VERBOSE("Solving the linear system ...\n");
        START_TIMING(time_solve);
        ExaGeoStatTrsmTile(EXAGEOSTAT_LEFT, EXAGEOSTAT_LOWER, EXAGEOSTAT_NO_TRANS, EXAGEOSTAT_NON_UNIT, 1, EXA_descriptorC, EXA_descriptorZ);
        STOP_TIMING(time_solve);
        flops = flops + FLOPS_DTRSM(ChamLeft, n, nhrs);
        VERBOSE(" Done.\n")

        //Calculate MLE likelihood
        VERBOSE("Calculating the MLE likelihood function ...");

        ExaGeoStatGemmTile(EXAGEOSTAT_TRANS, EXAGEOSTAT_NO_TRANS, 1, EXA_descriptorZ, EXA_descriptorZ, 0, EXA_descriptor_product);

        if (configurations->GetKernelName() == "bivariate_matern_parsimonious2_profile") {

            loglik = -(n / 2) + (n / 2) * log(n) - (n / 2) * log(*configurations->GetDotProduct()[0]) - 0.5 * logdet -
                     (T) (n / 2.0) * log(2.0 * PI);
            ExaGeoStatGemmTile(EXAGEOSTAT_TRANS, EXAGEOSTAT_NO_TRANS, 1, EXA_descZ1, EXA_descZ1, 0, EXA_descproduct1);
            ExaGeoStatGemmTile(EXAGEOSTAT_TRANS, EXAGEOSTAT_NO_TRANS, 1, EXA_descZ2, EXA_descZ2, 0, EXA_descproduct2);

            ///where do dotp1 and 2 get initialized???
            *configurations->GetVariance()[1] = (1.0 / (n / 2)) * (*configurations->GetDotProduct()[1]);
            *configurations->GetVariance()[2]  = (1.0 / (n / 2)) * (*configurations->GetDotProduct()[2]);

        }else {
            loglik = -0.5 * (*configurations->GetDotProduct()[1]) - 0.5 * logdet - (double) (n / 2.0) * log(2.0 * PI);
            *(configurations->GetVariance()[0]) = theta[0];
        }
        VERBOSE(" Done.\n");
    }

    fprintf(stderr, " %3d- Model Parameters (", iter_count + 1);
    if (configurations->GetLogger())
        ///should we be leaving the file opening and closing to the user?
        fprintf(configurations->GetFileLog(), " %3d- Model Parameters (",iter_count + 1);

    if ((configurations->GetKernelName() == "bivariate_matern_parsimonious_profile") ||
        (configurations->GetKernelName() == "bivariate_matern_parsimonious2_profile")) {
        fprintf(stderr, "%.8f, %.8f,", *configurations->GetVariance()[1], *configurations->GetVariance()[2]);
        if (configurations->GetLogger())
            fprintf(configurations->GetFileLog(), "%.8f, %.8f,", *configurations->GetVariance()[1], *configurations->GetVariance()[2]);
        i = 2;
    } else {
        i = 0;
    }
    for (; i < num_params; i++) {
        fprintf(stderr, "%.8f", theta[i]);
        if (i < num_params - 1) {
            fprintf(stderr, ",");
        }

        if (configurations->GetLogger()) {
            fprintf(configurations->GetFileLog(), "%.8f, ", theta[i]);
        }
    }
    fprintf(stderr, ")----> LogLi: %.18f\n", loglik);
    if (configurations->GetLogger()) {
        fprintf(configurations->GetFileLog(), ")----> LogLi: %.18f\n", loglik);
    }


    fprintf(stderr, " ---- Facto Time: %6.2f seconds\n", time_facto);
    fprintf(stderr, " ---- Matrix Generation Time: %6.2f\n", matrix_gen_time);
    fprintf(stderr, " ---- Total Time: %6.2f\n", matrix_gen_time + time_facto + logdet_calculate + time_solve);

    configurations->SetIterationsValue(configurations->GetIterationsValue() + 1);
    // for experiments
    if(Configurations::GetRunMode() == VERBOSE_MODE) {
        configurations->SetAvgExecutedTimePerIteration( configurations->GetAvgExecutedTimePerIteration() +/*matrix_gen_time*/+time_facto + logdet_calculate + time_solve);
        configurations->SetAvgFlopsPerIteration(configurations->GetAvgFlopsPerIteration() + flops / 1e9 / (time_facto + time_solve));
    }
     configurations->SetFinalLogLik(loglik);

    //output
////    results.final_loglik = loglik;
    return loglik;
}

template<typename T>
int ChameleonImplementationDense<T>::ExaGeoStatLapackCopyTile(exageostat::common::UpperLower aUpperLower, void *apA, void *apB){
    cout << "LPCY  descriptor m: " << ((CHAM_desc_t *) apA)->m << " n: " << ((CHAM_desc_t *) apA)->n << endl;
    cout << "LPCY B descriptor m: " << ((CHAM_desc_t *) apB)->m << " n: " << ((CHAM_desc_t *) apB)->n << endl;

    return CHAMELEON_dlacpy_Tile((cham_uplo_t) aUpperLower, (CHAM_desc_t *) apA, (CHAM_desc_t *) apB);
}

template<typename T>
int ChameleonImplementationDense<T>::ExaGeoStatLapackToDescriptor(exageostat::common::UpperLower aUpperLower, void *apAf77, int aLda, void * apA){
    return CHAMELEON_Lap2Desc((cham_uplo_t) aUpperLower, (CHAM_desc_t *) apAf77, aLda, (CHAM_desc_t *) apA);
}


template<typename T>
int ChameleonImplementationDense<T>::ExaGeoStatSequenceWait(void * apSequence) {
    return CHAMELEON_Sequence_Wait((RUNTIME_sequence_t *)apSequence);
}

template<typename T>
int ChameleonImplementationDense<T>::ExaGeoStatPotrfTile(exageostat::common::UpperLower aUpperLower, void *apA){
    return CHAMELEON_dpotrf_Tile((cham_uplo_t) aUpperLower, (CHAM_desc_t *) apA);
}

template<typename T>
int ChameleonImplementationDense<T>::ExaGeoStatTrsmTile(common::Side aSide, common::UpperLower aUpperLower, common::Trans aTrans, common::Diag aDiag, T aAlpha, void *apA, void *apB) {
    return CHAMELEON_dtrsm_Tile((cham_side_t) aSide, (cham_uplo_t) aUpperLower, (cham_trans_t) aTrans, (cham_diag_t) aDiag, aAlpha, (CHAM_desc_t *) apA, (CHAM_desc_t *) apB);
}

template<typename T>
int ChameleonImplementationDense<T>::ExaGeoStatGemmTile(common::Trans aTransA, common::Trans aTransB, T aAlpha,
                                                         void *apA, void *apB, T aBeta, void *apC) {
    return CHAMELEON_dgemm_Tile((cham_trans_t) aTransA, (cham_trans_t) aTransB, aAlpha,(CHAM_desc_t *) apA, (CHAM_desc_t *) apB, aBeta, (CHAM_desc_t *)apC);
}

template<typename T>
int ChameleonImplementationDense<T>::ExaGeoStatMeasureDetTileAsync(void *apDescA, void * apSequence, void *apRequest, void *apDescDet){
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
int ChameleonImplementationDense<T>::ExaGeoStaStrideVectorTileAsync(void *apDescA, void *apDescB, void *apDescC,
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
