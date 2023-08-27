
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ChameleonDescriptor.cpp
 * @brief Defines the ChameleonDescriptor class for creating matrix descriptors using the CHAMELEON library.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @author Mahmoud ElKarargy
 * @date 2023-08-15
**/

#include <data-units/descriptor/concrete/ChameleonDescriptor.hpp>

using namespace exageostat::dataunits::descriptor;

template<typename T>
CHAM_desc_t *ChameleonDescriptor<T>::CreateChameleonDescriptor(void *apDescriptor, const bool &aIsOOC, void *apMatrix,
                                                               const common::FloatPoint &aFloatPoint, const int &aMB,
                                                               const int &aNB, const int &aSize, const int &aLM,
                                                               const int &aLN, const int &aI, const int &aJ,
                                                               const int &aM, const int &aN, const int &aP,
                                                               const int &aQ) {
    auto chameleon_desc = (CHAM_desc_t *) apDescriptor;
    if (aIsOOC && apMatrix == nullptr && aMB != 1 && aNB != 1) {
        CHAMELEON_Desc_Create_OOC(&chameleon_desc, (cham_flttype_t) aFloatPoint, aMB, aNB, aSize, aLM, aLN, aI, aJ, aM,
                                  aN, aP, aQ);
    } else {
        CHAMELEON_Desc_Create(&chameleon_desc, apMatrix, (cham_flttype_t) aFloatPoint, aMB, aNB, aSize, aLM, aLN, aI,
                              aJ,
                              aM, aN, aP, aQ);
    }
    return chameleon_desc;
}

template<typename T>
int ChameleonDescriptor<T>::DestroyChameleonDescriptor(void *apDesc) {
    auto chameleon_desc = (CHAM_desc_t *) apDesc;
    return CHAMELEON_Desc_Destroy(&chameleon_desc);
}

