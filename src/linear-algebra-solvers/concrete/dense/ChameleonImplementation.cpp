
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

using namespace exageostat::linearAlgebra::dense;
using namespace exageostat::common;

template<typename T> void ChameleonImplementation<T>::InitiateDescriptors() {

    //// TODO: what C & Z stands for?
    CHAM_desc_t *pDescriptorC = nullptr;
    CHAM_desc_t *pDescriptorZ = nullptr;
    CHAM_desc_t *pDescriptorZcpy = nullptr;
    CHAM_desc_t *pDescriptorProduct = nullptr;
    CHAM_desc_t *pDescriptorDeterminant = nullptr;

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
    if (sizeof(T) == SIZE_OF_FLOAT){
        floatPoint = EXAGEOSTAT_REAL_FLOAT;
    }
    else {
        floatPoint = EXAGEOSTAT_REAL_DOUBLE;
    }


    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&pDescriptorC, nullptr, (cham_flttype_t) floatPoint, dts, dts, dts * dts, N, N, 0, 0, N, N, pGrid,
                                    qGrid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&pDescriptorZ, nullptr, (cham_flttype_t) floatPoint, dts, dts, dts * dts, N, 1, 0, 0, N, 1, pGrid,
                                    qGrid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&pDescriptorZcpy, Zcpy, (cham_flttype_t) floatPoint, dts, dts, dts * dts, N, 1, 0, 0, N, 1, pGrid,
                                    qGrid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&pDescriptorProduct, &dotProductValue, ChamRealFloat, dts, dts, dts * dts, 1, 1, 0, 0, 1,
                                    1, pGrid, qGrid);
    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&pDescriptorDeterminant, &dotProductValue, ChamRealFloat, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1,
                                    pGrid, qGrid);

}