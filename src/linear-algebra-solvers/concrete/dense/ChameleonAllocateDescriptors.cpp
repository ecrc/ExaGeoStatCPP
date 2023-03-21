
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

#include <linear-algebra-solvers/concrete/dense/ChameleonAllocateDescriptors.hpp>
#include <iostream>
#include <chameleon/struct.h>

using namespace exageostat::linearAlgebra::dense;
using namespace exageostat::dataunits;

template<typename T> void ChameleonAllocateDescriptors<T>::InitiateDescriptors(Precision aPrecision) {
    std::cout << "From Chameleon" << std::endl;

}

template<typename T> ChameleonAllocateDescriptors<T>::ChameleonAllocateDescriptors() {
    // Set the selected Precision
    std::cout << "Hi mahmoud, You selected: " << typeid(T).name() << " ." << std::endl;

}

template<typename T> void ChameleonAllocateDescriptors<T>::CreateDescriptors(T aPrecision) {

    //// TODO: what C & Z stands for?
    CHAM_desc_t *pDescriptorC = nullptr;
    CHAM_desc_t *pDescriptorZ = nullptr;
    CHAM_desc_t *pDescriptorZcpy = nullptr;
    CHAM_desc_t *pDescriptorProduct = nullptr;
    CHAM_desc_t *pDescriptorDeterminant = nullptr;

    RUNTIME_sequence_t *pSequence;
    RUNTIME_request_t request[2] = {CHAMELEON_SUCCESS, CHAMELEON_SUCCESS};

    // For ditributed system and should be removed
//    T *Zcpy = (float *) malloc(N * sizeof(float));
    T *Zcpy = (T *) malloc(1 * sizeof(T));

}
