
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file HicmaDescriptor.cpp
 * @brief Defines the Hicma Descriptor class for creating matrix descriptors using the HICMA library.
 * @version 1.0.1
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-08-15
**/
#include <data-units/descriptor/concrete/HicmaDescriptor.hpp>

using namespace exageostat::dataunits::descriptor;

template<typename T>
HICMA_desc_t *HicmaDescriptor<T>::CreateHicmaDescriptor(void *apDescriptor, const bool &aIsOOC, void *apMatrix,
                                                        const common::FloatPoint &aFloatPoint, const int &aMB,
                                                        const int &aNB, const int &aSize, const int &aLM,
                                                        const int &aLN, const int &aI, const int &aJ, const int &aM,
                                                        const int &aN, const int &aP, const int &aQ, const bool &aValidOOC) {
    auto hicma_desc = (HICMA_desc_t *) apDescriptor;
    if (aIsOOC && apMatrix == nullptr && aMB != 1 && aNB != 1 && aValidOOC) {
        HICMA_Desc_Create_OOC(&hicma_desc, (HICMA_enum) aFloatPoint, aMB, aNB, aSize, aLM, aLN, aI, aJ, aM, aN, aP, aQ);
    } else {
        HICMA_Desc_Create(&hicma_desc, apMatrix, (HICMA_enum) aFloatPoint, aMB, aNB, aSize, aLM, aLN, aI, aJ, aM, aN,
                          aP, aQ);
    }
    return hicma_desc;
}

template<typename T>
int HicmaDescriptor<T>::DestroyHicmaDescriptor(void *apDesc) {
    auto hicma_desc = (HICMA_desc_t *) apDesc;
    return HICMA_Desc_Destroy(&hicma_desc);
}
