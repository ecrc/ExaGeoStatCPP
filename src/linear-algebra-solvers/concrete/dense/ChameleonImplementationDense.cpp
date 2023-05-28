
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ChameleonAllocateDescriptors.cpp
 * @brief Sets up the Chameleon descriptors needed for the dense matrix computations in ExaGeoStat.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-03-20
**/

#include <linear-algebra-solvers/concrete/dense/ChameleonImplementationDense.hpp>
#include <lapacke.h>
// Include Chameleon libraries
extern "C" {
#include <chameleon/struct.h>
#include <chameleon.h>
#include <control/descriptor.h>
#include <control/context.h>
}

#define EXAGEOSTAT_RTBLKADDR(desc, type, m, n) ( (starpu_data_handle_t)RUNTIME_data_getaddr( desc, m, n ) )

// Use the following namespaces for convenience
using namespace exageostat::linearAlgebra::dense;
using namespace exageostat::common;
using namespace exageostat::dataunits;
using namespace exageostat::kernels;
using namespace std;

// Define a method to set up the Chameleon descriptors
template<typename T>
void ChameleonImplementationDense<T>::InitiateDescriptors() {

    // Initialize the Chameleon context
    this->ExaGeoStatInitContext(this->mpConfigurations->GetCoresNumber(), this->mpConfigurations->GetGPUsNumber());

    // Declare variables for Chameleon descriptors
    vector<void *> &pDescriptorC = this->mpConfigurations->GetDescriptorC();
    vector<void *> &pDescriptorZ = this->mpConfigurations->GetDescriptorZ();
    auto pDescriptorZcpy = (CHAM_desc_t **) &this->mpConfigurations->GetDescriptorZcpy();
    vector<void *> &pDescriptorProduct = this->mpConfigurations->GetDescriptorProduct();
    auto pDescriptorDeterminant = (CHAM_desc_t **) &this->mpConfigurations->GetDescriptorDeterminant();

    // Get the problem size and other configuration parameters
    int vectorSize;
    int N = this->mpConfigurations->GetProblemSize() * this->mpConfigurations->GetP();
    int dts = this->mpConfigurations->GetDenseTileSize();
    int pGrid = this->mpConfigurations->GetPGrid();
    int qGrid = this->mpConfigurations->GetQGrid();
    bool isOOC = this->mpConfigurations->GetIsOOC();

    // For distributed system and should be removed
    T *Zcpy = (T *) malloc(N * sizeof(T));
    T dotProductValue;

    // Create a Chameleon sequence
    RUNTIME_sequence_t *pSequence;
    RUNTIME_request_t request[2] = {CHAMELEON_SUCCESS, CHAMELEON_SUCCESS};
    CHAMELEON_Sequence_Create(&pSequence);

    // Set the floating point precision based on the template type
    FloatPoint floatPoint;
    if (sizeof(T) == SIZE_OF_FLOAT) {
        floatPoint = EXAGEOSTAT_REAL_FLOAT;
        vectorSize = 1;
    } else {
        floatPoint = EXAGEOSTAT_REAL_DOUBLE;
        vectorSize = 3;
    }

    // Create the Chameleon descriptors based on the configuration parameters
    // Depending on the passed Precession, the descriptor will resize for value 1 or 3.
    pDescriptorC.resize(vectorSize + 1, nullptr);
    pDescriptorZ.resize(vectorSize, nullptr);
    pDescriptorProduct.resize(vectorSize, nullptr);

    auto **CHAM_descriptorC = (CHAM_desc_t **) &pDescriptorC[0];
    EXAGEOSTAT_ALLOCATE_DENSE_MATRIX_TILE(CHAM_descriptorC, isOOC, nullptr, (cham_flttype_t) floatPoint, dts, dts,
                                          dts * dts, N, N, 0, 0, N, N, pGrid, qGrid)

    if (vectorSize > 1) {

        auto **CHAM_descsubC11 = (CHAM_desc_t **) &pDescriptorC[1];
        auto **CHAM_descsubC12 = (CHAM_desc_t **) &pDescriptorC[2];
        auto **CHAM_descsubC22 = (CHAM_desc_t **) &pDescriptorC[3];

        *CHAM_descsubC11 = chameleon_desc_submatrix(*CHAM_descriptorC, 0, 0, (*CHAM_descriptorC)->m / 2,
                                                    (*CHAM_descriptorC)->n / 2);
        *CHAM_descsubC12 = chameleon_desc_submatrix(*CHAM_descriptorC, (*CHAM_descriptorC)->m / 2, 0,
                                                    (*CHAM_descriptorC)->m / 2, (*CHAM_descriptorC)->n / 2);
        *CHAM_descsubC22 = chameleon_desc_submatrix(*CHAM_descriptorC, (*CHAM_descriptorC)->m / 2,
                                                    (*CHAM_descriptorC)->n / 2,
                                                    (*CHAM_descriptorC)->m / 2, (*CHAM_descriptorC)->n / 2);
    }
    EXAGEOSTAT_ALLOCATE_DENSE_MATRIX_TILE(((CHAM_desc_t **) &pDescriptorZ[0]), isOOC, nullptr,
                                          (cham_flttype_t) floatPoint, dts, dts, dts * dts, N, 1, 0, 0, N, 1, pGrid,
                                          qGrid)
    EXAGEOSTAT_ALLOCATE_DENSE_MATRIX_TILE(pDescriptorZcpy, isOOC, Zcpy, (cham_flttype_t) floatPoint, dts, dts,
                                          dts * dts, N, 1, 0, 0, N, 1, pGrid, qGrid)
    EXAGEOSTAT_ALLOCATE_DENSE_MATRIX_TILE(pDescriptorDeterminant, isOOC, &dotProductValue, (cham_flttype_t) floatPoint,
                                          dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, pGrid, qGrid)

    for (int idx = 1; idx < pDescriptorZ.size(); idx++) {
        auto **CHAM_descriptorZ_ = (CHAM_desc_t **) &pDescriptorZ[idx];
        EXAGEOSTAT_ALLOCATE_DENSE_MATRIX_TILE(CHAM_descriptorZ_, isOOC, nullptr, (cham_flttype_t) floatPoint, dts, dts,
                                              dts * dts, N / 2, 1, 0, 0, N / 2, 1, pGrid, qGrid)
    }

    for (auto &idx: pDescriptorProduct) {
        auto **CHAM_descriptorProduct = (CHAM_desc_t **) &idx;
        EXAGEOSTAT_ALLOCATE_DENSE_MATRIX_TILE(CHAM_descriptorProduct, isOOC, &dotProductValue,
                                              (cham_flttype_t) floatPoint, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, pGrid,
                                              qGrid)
    }

    this->mpConfigurations->SetSequence(pSequence);
    this->mpConfigurations->SetRequest(request);

    //stop gsl error handler
    gsl_set_error_handler_off();
}

template<typename T>
void ChameleonImplementationDense<T>::ExaGeoStatInitContext(const int &apCoresNumber, const int &apGPUs) {

    CHAM_context_t *chameleonContext;
    chameleonContext = chameleon_context_self();
    if (chameleonContext != nullptr) {
        printf("Another instance of Chameleon is already running...!");
    } else {
        CHAMELEON_user_tag_size(31, 26);
        CHAMELEON_Init(apCoresNumber, apGPUs)
    }
}

template<typename T>
void ChameleonImplementationDense<T>::ExaGeoStatFinalizeContext() {
    CHAM_context_t *chameleonContext;
    chameleonContext = chameleon_context_self();
    if (chameleonContext == nullptr) {
        printf("No active instance oh Chameleon...please use ExaGeoStatInitContext() function to initiate a new instance!\n");
    } else{
        CHAMELEON_Finalize();
    }
}

#define starpu_mpi_codelet(_codelet_) _codelet_

static void cl_dcmg_cpu_func(void *buffers[], void *cl_arg) {

    int m, n, m0, n0;
    exageostat::dataunits::Locations *apLocation1;
    exageostat::dataunits::Locations *apLocation2;
    exageostat::dataunits::Locations *apLocation3;
    double *theta;
    double *A;
    int distance_metric;
    exageostat::kernels::Kernel *kernel;

    A = (double *) STARPU_MATRIX_GET_PTR(buffers[0]);

    starpu_codelet_unpack_args(cl_arg, &m, &n, &m0, &n0, &apLocation1, &apLocation2, &apLocation3, &theta,
                               &distance_metric, &kernel);
    kernel->GenerateCovarianceMatrix(A, m, n, m0, n0, apLocation1,
                                     apLocation2, apLocation3, theta, distance_metric);
};


static struct starpu_codelet cl_dcmg =
        {
                .where        = STARPU_CPU /*| STARPU_CUDA*/,
                .cpu_func     = cl_dcmg_cpu_func,
#if defined(EXAGEOSTAT_USE_CUDA)
                //    .cuda_func      = {cl_dcmg_cuda_func},
#endif
                .nbuffers     = 1,
                .modes        = {STARPU_W},
                .name         = "dcmg"
        };

template<typename T>
void ChameleonImplementationDense<T>::CovarianceMatrixCodelet(void *descA, int uplo, dataunits::Locations *apLocation1,
                                                              dataunits::Locations *apLocation2,
                                                              dataunits::Locations *apLocation3,
                                                              double *aLocalTheta, int aDistanceMetric,
                                                              exageostat::kernels::Kernel *apKernel) {
    CHAM_context_t *chamctxt = chameleon_context_self();
    RUNTIME_option_t options;

    RUNTIME_options_init(&options, chamctxt, (RUNTIME_sequence_t *) this->mpConfigurations->GetSequence(),
                         (RUNTIME_request_t *) this->mpConfigurations->GetRequest());


    int tempmm, tempnn;

    auto *CHAM_descA = (CHAM_desc_t *) descA;
    CHAM_desc_t A = *CHAM_descA;
    struct starpu_codelet *cl = &cl_dcmg;
    int m, n, m0, n0;

    int size = A.n;

    for (n = 0; n < A.nt; n++) {
        tempnn = n == A.nt - 1 ? A.n - n * A.nb : A.nb;
        if (uplo == ChamUpperLower) {
            m = 0;
        } else {
            m = A.m == A.n ? n : 0;
        }
        for (; m < A.mt; m++) {

            tempmm = m == A.mt - 1 ? A.m - m * A.mb : A.mb;
            m0 = m * A.mb;
            n0 = n * A.nb;

            // Register the data with StarPU
            starpu_insert_task(starpu_mpi_codelet(cl),
                               STARPU_VALUE, &tempmm, sizeof(int),
                               STARPU_VALUE, &tempnn, sizeof(int),
                               STARPU_VALUE, &m0, sizeof(int),
                               STARPU_VALUE, &n0, sizeof(int),
                               STARPU_W, EXAGEOSTAT_RTBLKADDR(CHAM_descA, ChamRealDouble, m, n),
                               STARPU_VALUE, &apLocation1, sizeof(dataunits::Locations *),
                               STARPU_VALUE, &apLocation2, sizeof(dataunits::Locations *),
                               STARPU_VALUE, &apLocation3, sizeof(dataunits::Locations *),
                               STARPU_VALUE, &aLocalTheta, sizeof(double *),
                               STARPU_VALUE, &aDistanceMetric, sizeof(int),
                               STARPU_VALUE, &apKernel, sizeof(exageostat::kernels::Kernel *),
                               0);

            auto handle = EXAGEOSTAT_RTBLKADDR(CHAM_descA, ChamRealDouble, m, n);
            this->apMatrix = (double *) starpu_variable_get_local_ptr(handle);
        }
    }
    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, chamctxt);
}


template<typename T>
void ChameleonImplementationDense<T>::GenerateObservationsVector(void *descA, Locations *apLocation1,
                                                                 Locations *apLocation2, Locations *apLocation3,
                                                                 vector<double> aLocalTheta, int aDistanceMetric,
                                                                 Kernel *apKernel) {

    auto *sequence = (RUNTIME_sequence_t *) this->mpConfigurations->GetSequence();
    auto *request = (RUNTIME_request_t *) this->mpConfigurations->GetRequest();
    int N = this->mpConfigurations->GetProblemSize();

    //// TODO: Make all zeros, Seed.
    int iseed[4] = {0, 0, 0, 1};
    //nomral random generation of e -- ei~N(0, 1) to generate Z
    auto *Nrand = (double *) malloc(N * sizeof(double));
    LAPACKE_dlarnv(3, iseed, N, Nrand);

    //Generate the co-variance matrix C
//    VERBOSE("Initializing Covariance Matrix (Synthetic Dataset Generation Phase).....");
    auto *theta = (double *) malloc(aLocalTheta.size() * sizeof(double));
    for (int i = 0; i < aLocalTheta.size(); i++) {
        theta[i] = aLocalTheta[i];
    }
    this->CovarianceMatrixCodelet(descA, EXAGEOSTAT_LOWER, apLocation1, apLocation2, apLocation3, theta,
                                  aDistanceMetric, apKernel);

    CHAMELEON_Sequence_Wait(sequence);
    free(theta);
    //    VERBOSE(" Done.\n");

    //Copy Nrand to Z
//    VERBOSE("Generate Normal Random Distribution Vector Z (Synthetic Dataset Generation Phase) .....");
//    auto **CHAM_descriptorZ = (CHAM_desc_t **) &this->mpConfigurations->GetDescriptorZ()[0];
//    EXAGEOSTAT_Zcpy(*CHAM_descriptorZ, Nrand, sequence, request);
//    VERBOSE(" Done.\n");

    //Cholesky factorization for the Co-variance matrix C
//    VERBOSE("Cholesky factorization of Sigma (Synthetic Dataset Generation Phase) .....");
//    int success = CHAMELEON_dpotrf_Tile(ChamLower, (CHAM_desc_t *)descA);
//    SUCCESS(success, "Factorization cannot be performed..\n The matrix is not positive definite\n\n");
//    VERBOSE(" Done.\n");

    //Triangular matrix-matrix multiplication
//    VERBOSE("Triangular matrix-matrix multiplication Z=L.e (Synthetic Dataset Generation Phase) .....");
//    CHAMELEON_dtrmm_Tile(ChamLeft, ChamLower, ChamNoTrans, ChamNonUnit, 1, (CHAM_desc_t *) descA, *CHAM_descriptorZ);
//    VERBOSE(" Done.\n");

    //// TODO: make verbose in modes, Add log with path
//    if (log == 1) {
//        double *z;
//        CHAM_desc_t *CHAM_descZ = (CHAM_desc_t *) (data->descZ);
//        VERBOSE("Writing generated data to the disk (Synthetic Dataset Generation Phase) .....");
//#if defined(CHAMELEON_USE_MPI)
//        z = (double*) malloc(n * sizeof(double));
//        CHAMELEON_Tile_to_Lapack( CHAM_descZ, z, n);
//        if ( CHAMELEON_My_Mpi_Rank() == 0 )
//            write_vectors(z, data, n);
//        free(z);
//#else
//        z = CHAM_descZ->mat;
//        write_vectors(z, data, n);
    //free(z);
//#endif
//        VERBOSE(" Done.\n");
//    }
//
//    CHAMELEON_dlaset_Tile(ChamUpperLower, 0, 0, (CHAM_desc_t *) descA);
//    VERBOSE("Done Z Vector Generation Phase. (Chameleon Synchronous)\n");
//    VERBOSE("************************************************************\n");

}

template<typename T>
void
ChameleonImplementationDense<T>::EXAGEOSTAT_Zcpy(CHAM_desc_t *apDescA, double *apDoubleVector,
                                                 RUNTIME_sequence_t *apSequence,
                                                 RUNTIME_request_t *apRequest) {
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;

    chamctxt = chameleon_context_self();
    if (apSequence->status != CHAMELEON_SUCCESS)
        throw runtime_error("INVALID!");
    RUNTIME_options_init(&options, chamctxt, apSequence, apRequest);

    int m, m0;
    int tempmm;
    auto A = apDescA;
//    struct starpu_codelet *cl = &cl_dzcpy;

    for (m = 0; m < A->mt; m++) {
        tempmm = m == A->mt - 1 ? A->m - m * A->mb : A->mb;
        m0 = m * A->mb;

//        starpu_insert_task(starpu_mpi_codelet(cl),
//                           STARPU_VALUE, &tempmm, sizeof(int),
//                           STARPU_VALUE, &m0, sizeof(int),
//                           STARPU_VALUE, &apDoubleVector, sizeof(double),
//                           STARPU_W, RUNTIME_data_getaddr(((void *)A), m, 0),
//#if defined(CHAMELEON_CODELETS_HAVE_NAME)
//                STARPU_NAME, "dzcpy",
//#endif
//                           0);
//        core_dzcpy(((double *) RUNTIME_data_getaddr((A), m, 0)), tempmm, m0, apRequest);
        memcpy(((double *) RUNTIME_data_getaddr((A), m, 0)), &apRequest[m0], m * sizeof(double));
    }
    RUNTIME_options_ws_free(&options);
}