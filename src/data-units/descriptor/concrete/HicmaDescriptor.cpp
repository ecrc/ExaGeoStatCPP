/**
 * @file HicmaDescriptor.cpp
 * @brief Defines the Hicma Descriptor class for creating matrix descriptors using the HICMA library.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-08-15
**/
#include <data-units/descriptor/concrete/HicmaDescriptor.hpp>

using namespace exageostat::dataunits::descriptor;

template<typename T>
HICMA_desc_t *HicmaDescriptor<T>::CreateHicmaDescriptor(void *apDescriptor, bool aIsOOC, void *apMatrix,
                                                        common::FloatPoint aFloatPoint, int aMB, int aNB,
                                                        int aSize, int aLM, int aLN, int aI, int aJ, int aM,
                                                        int aN, int aP, int aQ) {
    auto hicma_desc = (HICMA_desc_t *) apDescriptor;
    if (aIsOOC && apMatrix == nullptr && aMB != 1 && aNB != 1) {
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
