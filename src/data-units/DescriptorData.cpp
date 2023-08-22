
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file DescriptorData.cpp
 * @brief Contains the definition of the DescriptorData class.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @author Mahmoud ElKarargy
 * @date 2023-07-18
**/

#include <data-units/DescriptorData.hpp>

using namespace exageostat::dataunits;
using namespace exageostat::common;
using namespace exageostat::dataunits::descriptor;

template<typename T>
DescriptorData<T>::DescriptorData(hardware::ExaGeoStatHardware &apHardware) {
    this->mpContext = apHardware.GetContext();
    if (!this->mpContext) {
        throw std::runtime_error("Can create descriptors, Hardware is not initialised!");
    }
}

template<typename T>
DescriptorData<T>::~DescriptorData() {

    ExaGeoStatDescriptor<T> exaGeoStatDescriptor;
    // Destroy descriptors.
    for (const auto &pair: this->mDictionary) {
        const std::string &key = pair.first;
        if (key.find("CHAMELEON") != std::string::npos && pair.second != nullptr) {
            exaGeoStatDescriptor.DestroyDescriptor(common::CHAMELEON_DESCRIPTOR, pair.second);
        } else if (pair.second != nullptr) {
            exaGeoStatDescriptor.DestroyDescriptor(common::HICMA_DESCRIPTOR, pair.second);
        }
    }
    this->mDictionary.clear();
#ifdef EXAGEOSTAT_USE_CHAMELEON
    if (this->mpSequence) {
        CHAMELEON_Sequence_Destroy((RUNTIME_sequence_t *) this->mpSequence);
    }
#endif
#ifdef EXAGEOSTAT_USE_HiCMA
    if(this->mpSequence){
        HICMA_Sequence_Destroy((HICMA_sequence_t *) this->mpSequence);
    }
#endif
}

template<typename T>
void *DescriptorData<T>::GetSequence() {
    return this->mpSequence;
}

template<typename T>
void DescriptorData<T>::SetRequest(void *apRequest) {
    this->mpRequest = apRequest;
}

template<typename T>
void DescriptorData<T>::SetSequence(void *apSequence) {
    this->mpSequence = apSequence;
}

template<typename T>
void *DescriptorData<T>::GetRequest() {
    return this->mpRequest;
}

template<typename T>
BaseDescriptor DescriptorData<T>::GetDescriptor(DescriptorType aDescriptorType, DescriptorName aDescriptorName) {

    BaseDescriptor descriptor{};
    if (aDescriptorType == CHAMELEON_DESCRIPTOR) {
#ifdef EXAGEOSTAT_USE_CHAMELEON
        if (this->mDictionary.find(GetDescriptorName(aDescriptorName) + "_CHAMELEON") == this->mDictionary.end()) {
            descriptor.chameleon_desc = nullptr;
        }
        descriptor.chameleon_desc = (CHAM_desc_t *) this->mDictionary[GetDescriptorName(aDescriptorName) +
                                                                      "_CHAMELEON"];
#else
        throw std::runtime_error("To use Chameleon descriptor you need to enable EXAGEOSTAT_USE_CHAMELEON!");
#endif

    } else {
#ifdef EXAGEOSTAT_USE_HiCMA
        descriptor.hicma_desc = (HICMA_desc_t *) this->mDictionary[GetDescriptorName(aDescriptorName) + "_HICMA"];
#else
        throw std::runtime_error("To use HiCMA descriptor you need to enable EXAGEOSTAT_USE_HiCMA!");
#endif
    }
    return descriptor;
}

template<typename T>
void DescriptorData<T>::SetDescriptor(DescriptorType aDescriptorType, DescriptorName aDescriptorName, bool aIsOOC,
                                      void *apMatrix, FloatPoint aFloatPoint, int aMB, int aNB, int aSize, int aLM,
                                      int aLN, int aI, int aJ, int aM, int aN, int aP, int aQ) {

    void *descriptor;
    std::string type;
    ExaGeoStatDescriptor<T> exaGeoStatDescriptor;
    if (aDescriptorType == CHAMELEON_DESCRIPTOR) {
#ifdef EXAGEOSTAT_USE_CHAMELEON
        descriptor = exaGeoStatDescriptor.CreateDescriptor((CHAM_desc_t *) descriptor, aDescriptorType, aIsOOC,
                                                           apMatrix, aFloatPoint, aMB, aNB, aSize, aLM, aLN, aI, aJ, aM,
                                                           aN, aP, aQ);
        type = "_CHAMELEON";
#else
        throw std::runtime_error("To create Chameleon descriptor you need to enable EXAGEOSTAT_USE_CHAMELEON!");
#endif
    } else {
#ifdef EXAGEOSTAT_USE_HiCMA
        descriptor = exaGeoStatDescriptor.CreateDescriptor((HICMA_desc_t *) descriptor, aDescriptorType, aIsOOC, apMatrix, aFloatPoint, aMB, aNB, aSize, aLM, aLN, aI, aJ, aM, aN, aP, aQ);
        type = "_HICMA";
#else
        throw std::runtime_error("To create HiCMA descriptor you need to enable EXAGEOSTAT_USE_HiCMA!");
#endif
    }

    this->mDictionary[GetDescriptorName(aDescriptorName) + type] = descriptor;
}

template<typename T>
T *DescriptorData<T>::GetDescriptorMatrix(common::DescriptorType aDescriptorType, void *apDesc) {
    if (aDescriptorType == common::CHAMELEON_DESCRIPTOR) {

#ifdef EXAGEOSTAT_USE_CHAMELEON
        return (T *) ((CHAM_desc_t *) apDesc)->mat;
#else
        throw std::runtime_error("To use Chameleon descriptor you need to enable EXAGEOSTAT_USE_CHAMELEON!");
#endif
    } else {
#ifdef EXAGEOSTAT_USE_HiCMA
        return (T *) ((HICMA_desc_t *) apDesc)->mat;
#else
        throw std::runtime_error("To use Hicma descriptor you need to enable EXAGEOSTAT_USE_HiCMA!");
#endif
    }
}

// Define a function that returns the name of a DescriptorName value as a string
template<typename T>
std::string DescriptorData<T>::GetDescriptorName(DescriptorName aDescriptorName) {
    switch (aDescriptorName) {
        case DESCRIPTOR_C:
            return "DESCRIPTOR_C";
        case DESCRIPTOR_C11:
            return "DESCRIPTOR_C11";
        case DESCRIPTOR_C12:
            return "DESCRIPTOR_C12";
        case DESCRIPTOR_C22:
            return "DESCRIPTOR_C22";
        case DESCRIPTOR_Z:
            return "DESCRIPTOR_Z";
        case DESCRIPTOR_Z_1:
            return "DESCRIPTOR_Z_1";
        case DESCRIPTOR_Z_2:
            return "DESCRIPTOR_Z_2";
        case DESCRIPTOR_Z_COPY:
            return "DESCRIPTOR_Z_COPY";
        case DESCRIPTOR_PRODUCT:
            return "DESCRIPTOR_PRODUCT";
        case DESCRIPTOR_PRODUCT_1:
            return "DESCRIPTOR_PRODUCT_1";
        case DESCRIPTOR_PRODUCT_2:
            return "DESCRIPTOR_PRODUCT_2";
        case DESCRIPTOR_DETERMINANT:
            return "DESCRIPTOR_DETERMINANT";
        case DESCRIPTOR_CD:
            return "DESCRIPTOR_CD";
        case DESCRIPTOR_CUV:
            return "DESCRIPTOR_CUV";
        case DESCRIPTOR_CRK:
            return "DESCRIPTOR_CRK";
        case DESCRIPTOR_Z_OBSERVATIONS:
            return "DESCRIPTOR_Z_OBSERVATIONS";
        case DESCRIPTOR_Z_Actual:
            return "DESCRIPTOR_Z_Actual";
        case DESCRIPTOR_MSE:
            return "DESCRIPTOR_MSE";
        case DESCRIPTOR_C12D:
            return "DESCRIPTOR_C12D";
        case DESCRIPTOR_C12UV:
            return "DESCRIPTOR_C12UV";
        case DESCRIPTOR_C12RK:
            return "DESCRIPTOR_C12RK";
        case DESCRIPTOR_C22D:
            return "DESCRIPTOR_C22D";
        case DESCRIPTOR_C22UV:
            return "DESCRIPTOR_C22UV";
        case DESCRIPTOR_C22RK:
            return "DESCRIPTOR_C22RK";
        default:
            throw std::invalid_argument(
                    "The name of descriptor you provided is undefined, Please read the user manual to know the available descriptors");
    }
}