
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

#include <linear-algebra-solvers/concrete/diagonal-super-tile/ChameleonImplementation.hpp>
#include <linear-algebra-solvers/concrete/dense/Helpers.hpp>
#include <chameleon/struct.h>
#include <chameleon.h>
#include <gsl/gsl_errno.h>

using namespace exageostat::linearAlgebra::diagonalSuperTile;
using namespace exageostat::common;
using namespace std;

template<typename T>
void ChameleonImplementation<T>::InitiateDescriptors() {

    vector<void *> pDescriptorC =  this->mpConfigurations->GetDescriptorC();
    vector<void *> pDescriptorZ = this->mpConfigurations->GetDescriptorZ();
    CHAM_desc_t * pChameleonDescriptorZcpy = (CHAM_desc_t*) this->mpConfigurations->GetDescriptorZcpy();
    vector<void *> pDescriptorProduct = this->mpConfigurations->GetDescriptorProduct();
    CHAM_desc_t * pChameleonDescriptorDeterminant = (CHAM_desc_t*) this->mpConfigurations->GetDescriptorDeterminant();

    pDescriptorC.push_back(nullptr);
    CHAM_desc_t * pChameleonDescriptorC = (CHAM_desc_t*) pDescriptorC[0];

    pDescriptorZ.push_back(nullptr);
    CHAM_desc_t * pChameleonDescriptorZ = (CHAM_desc_t*) pDescriptorZ[0];

    int vectorSize = 1;
    FloatPoint floatPoint = EXAGEOSTAT_REAL_FLOAT;
    if (sizeof(T) == SIZE_OF_FLOAT) {
        floatPoint = EXAGEOSTAT_REAL_FLOAT;
    } else {
        floatPoint = EXAGEOSTAT_REAL_DOUBLE;
        vectorSize = 3;
    }

    RUNTIME_sequence_t *pSequence;
    RUNTIME_request_t request[2] = {CHAMELEON_SUCCESS, CHAMELEON_SUCCESS};

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

    EXAGEOSTAT_ALLOCATE_DENSE_MATRIX_TILE(&pChameleonDescriptorC, isOOC, nullptr, (cham_flttype_t) floatPoint, dts, dts, dts * dts, N, N, 0, 0, N, N, pGrid, qGrid);
    EXAGEOSTAT_ALLOCATE_DENSE_MATRIX_TILE(&pChameleonDescriptorZ, isOOC, nullptr, (cham_flttype_t) floatPoint, dts, dts, dts * dts, N, 1, 0, 0, N, 1, pGrid, qGrid);
    EXAGEOSTAT_ALLOCATE_DENSE_MATRIX_TILE(&pChameleonDescriptorZcpy, isOOC, nullptr, (cham_flttype_t) floatPoint, dts, dts, dts * dts, N, 1, 0, 0, N, 1, pGrid, qGrid);

    for(int idx =0; idx <vectorSize; idx++){
        pDescriptorProduct.push_back(nullptr);
        CHAM_desc_t* pChameleonDescriptorProduct = (CHAM_desc_t*) pDescriptorProduct[idx];
        EXAGEOSTAT_ALLOCATE_DENSE_MATRIX_TILE(&pChameleonDescriptorProduct, isOOC, &dotProductValue, (cham_flttype_t) floatPoint, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, pGrid, qGrid);

    }
    EXAGEOSTAT_ALLOCATE_DENSE_MATRIX_TILE(&pChameleonDescriptorDeterminant, isOOC, &dotProductValue, (cham_flttype_t) floatPoint, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, pGrid, qGrid);

    //stop gsl error handler
    gsl_set_error_handler_off();
}