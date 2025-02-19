
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

#if DEFAULT_RUNTIME

#include <data-units/descriptor/concrete/ChameleonDescriptor.hpp>

#ifdef USE_HICMA
#include <data-units/descriptor/concrete/HicmaDescriptor.hpp>
#endif
#else
#include <data-units/descriptor/concrete/ParsecDescriptor.hpp>
#endif

using namespace exageostat::common;
using namespace exageostat::dataunits::descriptor;

template<typename T>
void *
ExaGeoStatDescriptor<T>::CreateDescriptor(void *apDescriptor, const DescriptorType &aDescriptorType,
                                          const bool &aIsOOC, void *apMatrix, const FloatPoint &aFloatPoint,
                                          const int &aMB, const int &aNB, const int &aSize, const int &aLM,
                                          const int &aLN, const int &aI, const int &aJ, const int &aM,
                                          const int &aN, const int &aP, const int &aQ, const bool &aValidOOC) {

#if DEFAULT_RUNTIME
    if (aDescriptorType == CHAMELEON_DESCRIPTOR) {
        return ChameleonDescriptor<T>::CreateChameleonDescriptor(apDescriptor, aIsOOC, apMatrix, aFloatPoint, aMB, aNB,
                                                                 aSize, aLM, aLN, aI, aJ, aM, aN, aP, aQ, aValidOOC);
    }
#ifdef USE_HICMA
    if (aDescriptorType == HICMA_DESCRIPTOR) {
        return HicmaDescriptor<T>::CreateHicmaDescriptor(apDescriptor, aIsOOC, apMatrix, aFloatPoint, aMB, aNB, aSize,
                                                         aLM, aLN, aI, aJ, aM, aN, aP, aQ, aValidOOC);
    }
#endif
#else
    if (aDescriptorType == PARSEC_DESCRIPTOR) {
        return ParsecDescriptor<T>::CreateParsecDescriptor(apDescriptor);
    }
#endif
    std::cerr << "Error, please select the correct descriptor type!" << std::endl;
    return nullptr;
}

template<typename T>
int ExaGeoStatDescriptor<T>::DestroyDescriptor(const DescriptorType &aDescriptorType, void *apDesc) {

#if DEFAULT_RUNTIME
    if (aDescriptorType == CHAMELEON_DESCRIPTOR) {
        return ChameleonDescriptor<T>::DestroyChameleonDescriptor(apDesc);
    }
#ifdef USE_HICMA
    if (aDescriptorType == HICMA_DESCRIPTOR) {
        return HicmaDescriptor<T>::DestroyHicmaDescriptor(apDesc);
    }
#endif
#else
    if (aDescriptorType == PARSEC_DESCRIPTOR) {
        return ParsecDescriptor<T>::DestroyParsecDescriptor(apDesc);
    }
#endif
    std::cerr << "Error, please select the correct descriptor type!" << std::endl;
    return -1;
}