
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

extern "C" {
#include <chameleon/struct.h>
#include <chameleon.h>
#include <control/context.h>
}
using namespace exageostat::linearAlgebra::diagonalSuperTile;
using namespace exageostat::common;
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

template<typename T>
void ChameleonImplementationDST<T>::CovarianceMatrixCodelet(void *descA, int uplo,
                                                            dataunits::Locations *apLocation1,
                                                            dataunits::Locations *apLocation2,
                                                            dataunits::Locations *apLocation3,
                                                            double *apLocalTheta, int aDistanceMetric,
                                                            exageostat::kernels::Kernel *apKernel) {

    int tempmm, tempnn;
    auto **A = (CHAM_desc_t **) &descA;
//    struct starpu_codelet *cl = &cl_dcmg;
    int m, n, m0, n0;

    int size = (*A)->n;

    std::cout << "nt: " << (*A)->nt << " mt: " << (*A)->mt << std::endl;
    for (n = 0; n < (*A)->nt; n++) {
        tempnn = n == (*A)->nt - 1 ? (*A)->n - n * (*A)->nb : (*A)->nb;
        if (uplo == ChamUpperLower)
            m = 0;
        else
            m = (*A)->m == (*A)->n ? n : 0;
        for (; m < (*A)->mt; m++) {

            tempmm = m == (*A)->mt - 1 ? (*A)->m - m * (*A)->mb : (*A)->mb;
            m0 = m * (*A)->mb;
            n0 = n * (*A)->nb;
            printf("m: %d n: %d uplo: %d\n", m, n, uplo);
            printf("EXAGEOSTAT_RTBLKADDR(descA, ChamRealDouble, m, n): %p\n", RUNTIME_data_getaddr((*A), m, n));

            apKernel->GenerateCovarianceMatrix(((double *) RUNTIME_data_getaddr((*A), m, n)), tempmm, tempnn, m0, n0,
                                               apLocation1, apLocation2, apLocation3, apLocalTheta, aDistanceMetric);
//            starpu_insert_task(starpu_mpi_codelet(cl),
//                               STARPU_VALUE, &tempmm, sizeof(int),
//                               STARPU_VALUE, &tempnn, sizeof(int),
//                               STARPU_VALUE, &m0, sizeof(int),
//                               STARPU_VALUE, &n0, sizeof(int),
//                               STARPU_W, EXAGEOSTAT_RTBLKADDR(descA, ChamRealDouble, m, n),
//                               STARPU_VALUE, &l1, sizeof(location * ),
//                               STARPU_VALUE, &l2, sizeof(location * ),
//                               STARPU_VALUE, &lm, sizeof(location * ),
//                               STARPU_VALUE, &theta, sizeof(double* ),
//                               STARPU_VALUE, &distance_metric, sizeof(int),
//                               STARPU_VALUE, &kernel, sizeof(int),
//                               0);
        }
    }
}
