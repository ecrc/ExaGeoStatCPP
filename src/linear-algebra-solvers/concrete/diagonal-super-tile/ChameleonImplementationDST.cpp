
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ChameleonAllocateDescriptors.cpp
 *
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