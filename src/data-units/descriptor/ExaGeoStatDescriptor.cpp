
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ExaGeoStatDescriptor.cpp
 * @brief Implementation of creating matrix descriptors used in CHAMELEON and HiCMA libraries.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @author Mahmoud ElKarargy
 * @date 2023-07-17
**/

#ifdef EXAGEOSTAT_USE_CHAMELEON

#include <data-units/descriptor/concrete/ChameleonDescriptor.hpp>

#endif
#ifdef EXAGEOSTAT_USE_HICMA
#include <data-units/descriptor/concrete/HicmaDescriptor.hpp>
#endif

#include <data-units/descriptor/ExaGeoStatDescriptor.hpp>

using namespace exageostat::common;
using namespace exageostat::dataunits::descriptor;

template<typename T>
void *ExaGeoStatDescriptor<T>::CreateDescriptor(void *apDescriptor, DescriptorType aDescriptorType, bool aIsOOC,
                                                void *apMatrix, FloatPoint aFloatPoint, int aMB, int aNB, int aSize,
                                                int aLM, int aLN, int aI, int aJ, int aM, int aN, int aP, int aQ) {

    if (aDescriptorType == CHAMELEON_DESCRIPTOR) {
#ifdef EXAGEOSTAT_USE_CHAMELEON
        return ChameleonDescriptor<T>::CreateChameleonDescriptor(apDescriptor, aIsOOC, apMatrix, aFloatPoint, aMB, aNB,
                                                                 aSize, aLM, aLN, aI, aJ, aM, aN, aP, aQ);
#endif
    } else if (aDescriptorType == HICMA_DESCRIPTOR) {
#ifdef EXAGEOSTAT_USE_HICMA
        return HicmaDescriptor<T>::CreateHicmaDescriptor(apDescriptor, aIsOOC, apMatrix, aFloatPoint, aMB, aNB, aSize, aLM, aLN, aI, aJ, aM, aN, aP, aQ);
#endif
    }
    std::cerr << "Error, please select the correct descriptor type!" << std::endl;
    return nullptr;
}

template<typename T>
int ExaGeoStatDescriptor<T>::DestroyDescriptor(DescriptorType aDescriptorType, void *apDesc) {
    if (aDescriptorType == CHAMELEON_DESCRIPTOR) {
#ifdef EXAGEOSTAT_USE_CHAMELEON
        return ChameleonDescriptor<T>::DestroyChameleonDescriptor(apDesc);
#endif
    } else if (aDescriptorType == HICMA_DESCRIPTOR) {
#ifdef EXAGEOSTAT_USE_HICMA
        return HicmaDescriptor<T>::DestroyHicmaDescriptor(apDesc);
#endif
    }
    std::cerr << "Error, please select the correct descriptor type!" << std::endl;
    return -1;
}