
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ExaGeoStatDescriptor.cpp
 * @brief Implementation of creating matrix descriptors used in CHAMELEON and HiCMA libraries.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-07-17
**/

#include <iostream>

#include <data-units/descriptor/ExaGeoStatDescriptor.hpp>
#include <data-units/descriptor/concrete/ChameleonDescriptor.hpp>

#ifdef USE_HICMA

#include <data-units/descriptor/concrete/HicmaDescriptor.hpp>

#endif

using namespace exageostat::common;
using namespace exageostat::dataunits::descriptor;

template<typename T>
void *
ExaGeoStatDescriptor<T>::CreateDescriptor(void *apDescriptor, const common::DescriptorType &aDescriptorType,
                                          const bool &aIsOOC, void *apMatrix, const common::FloatPoint &aFloatPoint,
                                          const int &aMB, const int &aNB, const int &aSize, const int &aLM,
                                          const int &aLN, const int &aI, const int &aJ, const int &aM,
                                          const int &aN, const int &aP, const int &aQ, const bool &aValidOOC) {

    if (aDescriptorType == CHAMELEON_DESCRIPTOR) {
        return ChameleonDescriptor<T>::CreateChameleonDescriptor(apDescriptor, aIsOOC, apMatrix, aFloatPoint, aMB, aNB,
                                                                 aSize, aLM, aLN, aI, aJ, aM, aN, aP, aQ, aValidOOC);
    } else if (aDescriptorType == HICMA_DESCRIPTOR) {
#ifdef USE_HICMA
        return HicmaDescriptor<T>::CreateHicmaDescriptor(apDescriptor, aIsOOC, apMatrix, aFloatPoint, aMB, aNB, aSize,
                                                         aLM, aLN, aI, aJ, aM, aN, aP, aQ, aValidOOC);
#endif
    }
    std::cerr << "Error, please select the correct descriptor type!" << std::endl;
    return nullptr;
}

template<typename T>
int ExaGeoStatDescriptor<T>::DestroyDescriptor(const DescriptorType &aDescriptorType, void *apDesc) {
    if (aDescriptorType == CHAMELEON_DESCRIPTOR) {
        return ChameleonDescriptor<T>::DestroyChameleonDescriptor(apDesc);
    } else if (aDescriptorType == HICMA_DESCRIPTOR) {
#ifdef USE_HICMA
        return HicmaDescriptor<T>::DestroyHicmaDescriptor(apDesc);
#endif
    }
    std::cerr << "Error, please select the correct descriptor type!" << std::endl;
    return -1;
}