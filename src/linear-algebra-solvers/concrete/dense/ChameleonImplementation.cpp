
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
#include <iostream>
#include <chameleon/struct.h>

using namespace exageostat::linearAlgebra::dense;
using namespace exageostat::common;

template<typename T> void ChameleonImplementation<T>::InitiateDescriptors(T aPrecision) {

    //// TODO: what C & Z stands for?
    CHAM_desc_t *pDescriptorC = nullptr;
    CHAM_desc_t *pDescriptorZ = nullptr;
    CHAM_desc_t *pDescriptorZcpy = nullptr;
    CHAM_desc_t *pDescriptorProduct = nullptr;
    CHAM_desc_t *pDescriptorDeterminant = nullptr;

    RUNTIME_sequence_t *pSequence;
    RUNTIME_request_t request[2] = {CHAMELEON_SUCCESS, CHAMELEON_SUCCESS};

    // For ditributed system and should be removed
//    T *Zcpy = (T *) malloc(aProblemSize * sizeof(T));


}