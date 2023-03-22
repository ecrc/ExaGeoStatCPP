
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

#include <linear-algebra-solvers/concrete/dense/ChameleonImplementation.hpp>
#include <linear-algebra-solvers/concrete/dense/Helpers.hpp>
#include <iostream>
#include <chameleon/struct.h>
#include <chameleon.h>
#include <vector>

using namespace exageostat::linearAlgebra::dense;
using namespace exageostat::common;

template<typename T>
void ChameleonImplementation<T>::InitiateDescriptors() {

    //// TODO: what C & Z stands for?
    std::vector<CHAM_desc_t *> pDescriptorC;
    std::vector<CHAM_desc_t *> pDescriptorZ;
    CHAM_desc_t  *pDescriptorZcpy = nullptr;
    std::vector<CHAM_desc_t *> pDescriptorProduct;
    CHAM_desc_t *pDescriptorDeterminant = nullptr;
    int vectorSize = 0;
    
    RUNTIME_sequence_t *pSequence;
    RUNTIME_request_t request[2] = {CHAMELEON_SUCCESS, CHAMELEON_SUCCESS};

    int N = this->mpConfigurations->GetProblemSize() * this->mpConfigurations->GetP();
    int dts = this->mpConfigurations->GetTileSize();
    int pGrid = this->mpConfigurations->GetPGrid();
    int qGrid = this->mpConfigurations->GetQGrid();

    // For distributed system and should be removed
    T *Zcpy = (T *) malloc(N * sizeof(T));

    T dotProductValue;

    //Identifies a set of routines sharing common exception handling.
    CHAMELEON_Sequence_Create(&pSequence);

    FloatPoint floatPoint;
    if (sizeof(T) == SIZE_OF_FLOAT) {
        floatPoint = EXAGEOSTAT_REAL_FLOAT;
        vectorSize = 1;
    } else {
        floatPoint = EXAGEOSTAT_REAL_DOUBLE;
        vectorSize = 3;
    }
    // Depending on the passed Precession, A for loop with value 1 or 3 will initialize the vectors
    for (int idx = 0; idx < vectorSize; idx++){
        pDescriptorC.push_back(nullptr);
        pDescriptorZ.push_back(nullptr);
        pDescriptorProduct.push_back(nullptr);
    }

    if(vectorSize > 1){
        pDescriptorC.push_back(nullptr);
//        pDescriptorC[1] = chameleon_desc_submatrix(pDescriptorC[0], 0, 0, pDescriptorC[0]->m / 2, pDescriptorC[0]->n / 2);
//        pDescriptorC[2] = chameleon_desc_submatrix(pDescriptorC[0], pDescriptorC[0]->m / 2, 0, pDescriptorC[0]->m / 2, pDescriptorC[0]->n / 2);
//        pDescriptorC[3] = chameleon_desc_submatrix(pDescriptorC[0], pDescriptorC[0]->m / 2, pDescriptorC[0]->n / 2, pDescriptorC[0]->m / 2, pDescriptorC[0]->n / 2);

    }

    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&pDescriptorC[0], nullptr, (cham_flttype_t) floatPoint, dts, dts, dts * dts, N, N, 0, 0, N, N, pGrid, qGrid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&pDescriptorZ[0], nullptr, (cham_flttype_t) floatPoint, dts, dts, dts * dts, N, 1, 0,0, N, 1, pGrid, qGrid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&pDescriptorZcpy, Zcpy, (cham_flttype_t) floatPoint, dts, dts, dts * dts, N, 1, 0, 0, N, 1, pGrid, qGrid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&pDescriptorDeterminant, &dotProductValue, (cham_flttype_t) floatPoint, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, pGrid, qGrid);

    for (int idx = 1; idx < pDescriptorZ.size(); idx++) {
        EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&pDescriptorZ[idx], nullptr, (cham_flttype_t) floatPoint, dts, dts, dts * dts, N / 2, 1, 0, 0, N / 2, 1, pGrid, qGrid);
    }

    for (int idx = 0; idx < pDescriptorZ.size(); idx++) {
        EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&pDescriptorProduct[idx], &dotProductValue, (cham_flttype_t) floatPoint, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1, pGrid, qGrid)
    }

    //stop gsl error handler
//    gsl_set_error_handler_off();
}