
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

    int problemSize = this->mpConfigurations->GetProblemSize();

    // For distributed system and should be removed
    T *Zcpy = (T *) malloc(problemSize * sizeof(T));

    //Identifies a set of routines sharing common exception handling.
//    CHAMELEON_Sequence_Create(&pSequence);

    FloatPoint floatPoint;
    if (sizeof(T) == SIZE_OF_FLOAT){
        floatPoint = EXAGEOSTAT_REAL_FLOAT;
    }
    else {
        floatPoint = EXAGEOSTAT_REAL_DOUBLE;
    }

//    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&pDescriptorC, nullptr, floatPoint, dts, dts, dts * dts, N, N, 0, 0, N, N, p_grid,
//                                    q_grid);
//    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descZ, nullptr, floatPoint, dts, dts, dts * dts, N, 1, 0, 0, N, 1, p_grid,
//                                    q_grid);
//    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descZcpy, Zcpy, floatPoint, dts, dts, dts * dts, N, 1, 0, 0, N, 1, p_grid,
//                                    q_grid);
//    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descproduct, &data->sdotp, floatPoint, dts, dts, dts * dts, 1, 1, 0, 0, 1,
//                                    1, p_grid, q_grid);
//    EXAGEOSTAT_ALLOCATE_MATRIX_TILE(&CHAM_descdet, &data->det, floatPoint, dts, dts, dts * dts, 1, 1, 0, 0, 1, 1,
//                                    p_grid, q_grid);


}