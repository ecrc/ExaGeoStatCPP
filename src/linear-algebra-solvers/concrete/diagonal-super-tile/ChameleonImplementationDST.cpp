
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ChameleonAllocateDescriptors.cpp
 * @brief Sets up the Chameleon descriptors needed for the diagonal super tile computations in ExaGeoStat.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-03-20
**/

#include <linear-algebra-solvers/concrete/diagonal-super-tile/ChameleonImplementationDST.hpp>
#include <lapacke.h>

extern "C" {
#include <chameleon/struct.h>
#include <chameleon.h>
#include <control/context.h>
}
using namespace exageostat::linearAlgebra::diagonalSuperTile;
using namespace exageostat::common;
using namespace exageostat::kernels;
using namespace exageostat::dataunits;
using namespace std;

template<typename T>
void ChameleonImplementationDST<T>::InitiateDescriptors() {

    // Initialize Exageostat Hardware.
    this->ExaGeoStatInitContext(this->mpConfigurations->GetCoresNumber(), this->mpConfigurations->GetGPUsNumber());

    vector<void *> &pDescriptorC = this->mpConfigurations->GetDescriptorC();
    vector<void *> &pDescriptorZ = this->mpConfigurations->GetDescriptorZ();
    auto pChameleonDescriptorZcpy = (CHAM_desc_t **) &this->mpConfigurations->GetDescriptorZcpy();
    vector<void *> &pDescriptorProduct = this->mpConfigurations->GetDescriptorProduct();
    auto pChameleonDescriptorDeterminant = (CHAM_desc_t **) &this->mpConfigurations->GetDescriptorDeterminant();

    pDescriptorC.push_back(nullptr);
    auto **pChameleonDescriptorC = (CHAM_desc_t **) &pDescriptorC[0];

    pDescriptorZ.push_back(nullptr);
    auto **pChameleonDescriptorZ = (CHAM_desc_t **) &pDescriptorZ[0];

    int vectorSize = 1;
    FloatPoint floatPoint;
    if (sizeof(T) == SIZE_OF_FLOAT) {
        floatPoint = EXAGEOSTAT_REAL_FLOAT;
    } else {
        floatPoint = EXAGEOSTAT_REAL_DOUBLE;
        vectorSize = 3;
    }

    RUNTIME_sequence_t *pSequence;

    int N = this->mpConfigurations->GetProblemSize() * this->mpConfigurations->GetP();
    int dts = this->mpConfigurations->GetDenseTileSize();
    int pGrid = this->mpConfigurations->GetPGrid();
    int qGrid = this->mpConfigurations->GetQGrid();
    bool isOOC = this->mpConfigurations->GetIsOOC();

    // For distributed system and should be removed
    T *Zcpy = (T *) malloc(N * sizeof(T));
    T dotProductValue;

    //Identifies a set of routines sharing common exception handling.
    CHAMELEON_Sequence_Create(&pSequence);

    EXAGEOSTAT_ALLOCATE_DENSE_MATRIX_TILE(pChameleonDescriptorC, isOOC, nullptr, (cham_flttype_t) floatPoint, dts, dts,
                                          dts * dts, N, N, 0, 0, N, N, pGrid, qGrid)
    EXAGEOSTAT_ALLOCATE_DENSE_MATRIX_TILE(pChameleonDescriptorZ, isOOC, nullptr, (cham_flttype_t) floatPoint, dts, dts,
                                          dts * dts, N, 1, 0, 0, N, 1, pGrid, qGrid)
    EXAGEOSTAT_ALLOCATE_DENSE_MATRIX_TILE(pChameleonDescriptorZcpy, isOOC, Zcpy, (cham_flttype_t) floatPoint, dts,
                                          dts, dts * dts, N, 1, 0, 0, N, 1, pGrid, qGrid)

    for (int idx = 0; idx < vectorSize; idx++) {
        pDescriptorProduct.push_back(nullptr);
        auto **pChameleonDescriptorProduct = (CHAM_desc_t **) &pDescriptorProduct[idx];
        EXAGEOSTAT_ALLOCATE_DENSE_MATRIX_TILE(pChameleonDescriptorProduct, isOOC, &dotProductValue,
                                              (cham_flttype_t) floatPoint, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, pGrid,
                                              qGrid)

    }
    EXAGEOSTAT_ALLOCATE_DENSE_MATRIX_TILE(pChameleonDescriptorDeterminant, isOOC, &dotProductValue,
                                          (cham_flttype_t) floatPoint, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, pGrid,
                                          qGrid)
    this->ExaGeoStatFinalizeContext();
    //stop gsl error handler
    gsl_set_error_handler_off();
}

template<typename T>
void ChameleonImplementationDST<T>::ExaGeoStatInitContext(const int &apCoresNumber, const int &apGPUs) {

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
void ChameleonImplementationDST<T>::ExaGeoStatFinalizeContext() {

    CHAM_context_t *chameleonContext;
    chameleonContext = chameleon_context_self();
    if (chameleonContext == nullptr) {
        printf("No active instance oh Chameleon...please use ExaGeoStatInitContext() function to initiate a new instance!\n");
    } else
        CHAMELEON_Finalize();
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
#define EXAGEOSTAT_RTBLKADDR(desc, type, m, n) ( (starpu_data_handle_t)RUNTIME_data_getaddr( desc, m, n ) )

template<typename T>
void ChameleonImplementationDST<T>::CovarianceMatrixCodelet(void *descA, int uplo,
                                                            dataunits::Locations *apLocation1,
                                                            dataunits::Locations *apLocation2,
                                                            dataunits::Locations *apLocation3,
                                                            double *aLocalTheta, int aDistanceMetric,
                                                            exageostat::kernels::Kernel *apKernel) {
    CHAM_context_t *chamctxt;
    RUNTIME_option_t options;
    chamctxt = chameleon_context_self();

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
void ChameleonImplementationDST<T>::GenerateObservationsVector(void *descA, Locations *apLocation1,
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
}