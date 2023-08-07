
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ExaGeoStatDescriptor.cpp
 * @brief 
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-07-17
**/

#include <iostream>
#include <string>

#ifdef EXAGEOSTAT_USE_CHAMELEON
extern "C" {
#include <chameleon/struct.h>
#include <chameleon.h>
#include <control/descriptor.h>
}
#endif

#ifdef EXAGEOSTAT_USE_HiCMA
extern "C"{
#include <hicma.h>
}
#endif

#include <data-units/ExaGeoStatDescriptor.hpp>

using namespace exageostat::common;
using namespace exageostat::dataunits;

#ifdef EXAGEOSTAT_USE_CHAMELEON

template<typename T>
CHAM_desc_t *ExaGeoStatDescriptor<T>::CreateChameleonDescriptor(CHAM_desc_t *apDescriptor, bool aIsOOC, void *apMatrix,
                                                                FloatPoint aFloatPoint, int aMB, int aNB, int aSize,
                                                                int aLM, int aLN, int aI, int aJ, int aM, int aN,
                                                                int aP, int aQ) {

    if (aIsOOC && apMatrix == nullptr && aMB != 1 && aNB != 1) {
        CHAMELEON_Desc_Create_OOC(&apDescriptor, (cham_flttype_t) aFloatPoint, aMB, aNB, aSize, aLM, aLN, aI, aJ, aM,
                                  aN, aP, aQ);
    } else {
        CHAMELEON_Desc_Create(&apDescriptor, apMatrix, (cham_flttype_t) aFloatPoint, aMB, aNB, aSize, aLM, aLN, aI, aJ,
                              aM, aN, aP, aQ);
    }
    return apDescriptor;
}

template<typename T> CHAM_desc_t * ExaGeoStatDescriptor<T>::CreateChameleonSubMatrixDescriptor(CHAM_desc_t *apDescriptor, int aI, int aJ, int aM, int aN) {

    return chameleon_desc_submatrix(apDescriptor, aI, aJ, aM, aN);
}

#endif

#ifdef EXAGEOSTAT_USE_HiCMA

template<typename T>
HICMA_desc_t *ExaGeoStatDescriptor<T>::CreateHicmaDescriptor(HICMA_desc_t *apDescriptor, bool aIsOOC, void *apMatrix,
                                                             FloatPoint aFloatPoint, int aMB, int aNB, int aSize,
                                                             int aLM, int aLN, int aI, int aJ, int aM, int aN,
                                                             int aP, int aQ) {

    if (aIsOOC && apMatrix == nullptr && aMB != 1 && aNB != 1) {
        HICMA_Desc_Create_OOC(&apDescriptor, (HICMA_enum) aFloatPoint, aMB, aNB, aSize, aLM, aLN, aI, aJ, aM, aN, aP, aQ);
    } else {
        HICMA_Desc_Create(&apDescriptor, apMatrix, (HICMA_enum) aFloatPoint, aMB, aNB, aSize, aLM, aLN, aI, aJ, aM, aN, aP, aQ);
    }

    return apDescriptor;
}

#endif