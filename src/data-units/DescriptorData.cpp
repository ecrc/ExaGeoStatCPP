
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file DescriptorData.cpp
 * @brief Contains the definition of the DescriptorData class.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-07-18
**/

#include <data-units/DescriptorData.hpp>

using namespace exageostat::dataunits;
using namespace exageostat::common;
using namespace exageostat::dataunits::descriptor;

template<typename T>
DescriptorData<T>::~DescriptorData() {

    ExaGeoStatDescriptor<T> exageostat_descriptor;
    // Destroy descriptors.
    const std::string &chameleon = "_CHAMELEON";
    for (const auto &pair: this->mDictionary) {
        const std::string &key = pair.first;
#if DEFAULT_RUNTIME
        if (key.find("CHAMELEON") != std::string::npos && pair.second != nullptr) {
            exageostat_descriptor.DestroyDescriptor(CHAMELEON_DESCRIPTOR, pair.second);
#ifdef USE_HICMA
            // Since there are converted descriptors from Chameleon to Hicma, which have the same memory address.
            // So, by deleting the owner which is Chameleon, no need to delete hicma. Therefore, we remove the row of that descriptor.
            std::string converted_chameleon = key.substr(0, key.length() - chameleon.length());
            const std::string &desc = converted_chameleon + "_CHAM_HIC";
            if (this->mDictionary.find(desc) != this->mDictionary.end()) {
                delete (HICMA_desc_t *)
                        this->mDictionary[converted_chameleon + "_CHAM_HIC"];
                this->mDictionary.erase(converted_chameleon + "_CHAM_HIC");
            }
#endif
        } else if (key.find("HICMA") != std::string::npos && pair.second != nullptr) {
            exageostat_descriptor.DestroyDescriptor(HICMA_DESCRIPTOR, pair.second);
        }
#else
    if (key.find("PARSEC") != std::string::npos && pair.second != nullptr) {
            if(key == "DESCRIPTOR_FLMT_PARSEC") {
                continue;
            }
            exageostat_descriptor.DestroyDescriptor(PARSEC_DESCRIPTOR, pair.second);
    }
#endif

    }
    this->mDictionary.clear();
#if DEFAULT_RUNTIME
    if (this->mpSequence) {
        CHAMELEON_Sequence_Destroy((RUNTIME_sequence_t *) this->mpSequence);
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

#ifdef USE_HICMA

template<typename T>
HICMA_desc_t *
DescriptorData<T>::ConvertChameleonToHicma(CHAM_desc_t *apChameleonDesc, const DescriptorName &aDescriptorName) {

    this->SetDescriptor(HICMA_DESCRIPTOR, aDescriptorName, apChameleonDesc->ooc, apChameleonDesc->mat,
                        (FloatPoint) (apChameleonDesc->dtyp), apChameleonDesc->mb, apChameleonDesc->nb,
                        apChameleonDesc->bsiz, apChameleonDesc->lm, apChameleonDesc->ln, apChameleonDesc->i,
                        apChameleonDesc->j, apChameleonDesc->m, apChameleonDesc->n, apChameleonDesc->p,
                        apChameleonDesc->q, apChameleonDesc->ooc, true);

    return this->GetDescriptor(HICMA_DESCRIPTOR, aDescriptorName).hicma_desc;
}

#endif

template<typename T>
BaseDescriptor
DescriptorData<T>::GetDescriptor(const DescriptorType &aDescriptorType, const DescriptorName &aDescriptorName) {

    BaseDescriptor descriptor{};
#if DEFAULT_RUNTIME

    if (aDescriptorType == CHAMELEON_DESCRIPTOR) {
        if (this->mDictionary.find(GetDescriptorName(aDescriptorName) + "_CHAMELEON") == this->mDictionary.end()) {
            descriptor.chameleon_desc = nullptr;
        }
        descriptor.chameleon_desc = (CHAM_desc_t *) this->mDictionary[GetDescriptorName(aDescriptorName) +
                                                                      "_CHAMELEON"];
    } else {
#ifdef USE_HICMA
        if (this->mDictionary.find(GetDescriptorName(aDescriptorName) + "_HICMA") != this->mDictionary.end()) {
            descriptor.hicma_desc = (HICMA_desc_t *)
                    this->mDictionary[GetDescriptorName(aDescriptorName) + "_HICMA"];
        } else if (this->mDictionary.find(GetDescriptorName(aDescriptorName) + "_CHAM_HIC") !=
                   this->mDictionary.end()) {
            descriptor.hicma_desc = (HICMA_desc_t *)
                    this->mDictionary[GetDescriptorName(aDescriptorName) + "_CHAM_HIC"];
        } else if (this->mDictionary.find(GetDescriptorName(aDescriptorName) + "_CHAMELEON") !=
                   this->mDictionary.end()) {
            descriptor.hicma_desc = this->ConvertChameleonToHicma(
                    (CHAM_desc_t *) this->mDictionary[GetDescriptorName(aDescriptorName) + "_CHAMELEON"],
                    aDescriptorName);
            this->mDictionary[GetDescriptorName(aDescriptorName) + "_CHAM_HIC"] = descriptor.hicma_desc;
        } else {
            descriptor.hicma_desc = nullptr;
        }
#else
        throw std::runtime_error("To use HiCMA descriptor you need to enable USE_HICMA!");
#endif
    }
#else
    if (aDescriptorType == PARSEC_DESCRIPTOR) {
        if (this->mDictionary.find(GetDescriptorName(aDescriptorName) + "_PARSEC") == this->mDictionary.end()) {
            descriptor.parsec_desc = nullptr;
        }
        descriptor.parsec_desc = (parsec_matrix_block_cyclic_t *) this->mDictionary[GetDescriptorName(aDescriptorName) +
                                                                      "_PARSEC"];
    }
#endif
    return descriptor;
}

template<typename T>
void DescriptorData<T>::SetDescriptor(const DescriptorType &aDescriptorType, const DescriptorName &aDescriptorName,
                                      const bool &aIsOOC, void *apMatrix, const FloatPoint &aFloatPoint,
                                      const int &aMB, const int &aNB, const int &aSize, const int &aLM, const int &aLN,
                                      const int &aI, const int &aJ, const int &aM, const int &aN, const int &aP,
                                      const int &aQ, const bool &aValidOOC, const bool &aConverted) {

    void *descriptor;
    std::string type;
    ExaGeoStatDescriptor<T> exageostat_descriptor;
#if DEFAULT_RUNTIME

    if (aDescriptorType == CHAMELEON_DESCRIPTOR) {
        descriptor = exageostat_descriptor.CreateDescriptor((CHAM_desc_t *) descriptor, aDescriptorType, aIsOOC,
                                                           apMatrix, aFloatPoint, aMB, aNB, aSize, aLM, aLN, aI, aJ, aM,
                                                           aN, aP, aQ, aValidOOC);
        type = "_CHAMELEON";

    } else {
#ifdef USE_HICMA
        descriptor = exageostat_descriptor.CreateDescriptor((HICMA_desc_t *) descriptor, aDescriptorType, aIsOOC,
                                                           apMatrix, aFloatPoint, aMB, aNB, aSize, aLM, aLN, aI, aJ, aM,
                                                           aN, aP, aQ, aValidOOC);
        type = "_HICMA";
#else
        throw std::runtime_error("To create HiCMA descriptor you need to enable USE_HICMA!");
#endif
    }
    if (aConverted) {
        type = "_CHAM_HIC";
    }
#else
    if (aDescriptorType == PARSEC_DESCRIPTOR) {
        descriptor = exageostat_descriptor.CreateDescriptor((parsec_matrix_block_cyclic_t *) descriptor, aDescriptorType, aIsOOC,
                                                           apMatrix, aFloatPoint, aMB, aNB, aSize, aLM, aLN, aI, aJ, aM,
                                                           aN, aP, aQ, aValidOOC);
        type = "_PARSEC";
    }
    else {
        throw std::runtime_error("While using PaRSEC as a runtime, only PaRSEC descriptors are enabled!");
    }
#endif
    this->mDictionary[GetDescriptorName(aDescriptorName) + type] = descriptor;
}

template<typename T>
T *
DescriptorData<T>::GetDescriptorMatrix(const DescriptorType &aDescriptorType, const DescriptorName &aDescriptorName) {
#if DEFAULT_RUNTIME
    if (aDescriptorType == CHAMELEON_DESCRIPTOR) {
        return (T *) (this->GetDescriptor(CHAMELEON_DESCRIPTOR, aDescriptorName).chameleon_desc)->mat;
    } else {
#ifdef USE_HICMA
        return (T *) (this->GetDescriptor(HICMA_DESCRIPTOR, aDescriptorName).hicma_desc)->mat;
#else
        throw std::runtime_error("To use Hicma descriptor you need to enable USE_HICMA!");
#endif
    }
#endif
}

// Define a function that returns the name of a DescriptorName value as a string
template<typename T>
std::string DescriptorData<T>::GetDescriptorName(const DescriptorName &aDescriptorName) {
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
        case DESCRIPTOR_Z_3:
            return "DESCRIPTOR_Z_3";
        case DESCRIPTOR_Z_COPY:
            return "DESCRIPTOR_Z_COPY";
        case DESCRIPTOR_PRODUCT:
            return "DESCRIPTOR_PRODUCT";
        case DESCRIPTOR_PRODUCT_1:
            return "DESCRIPTOR_PRODUCT_1";
        case DESCRIPTOR_PRODUCT_2:
            return "DESCRIPTOR_PRODUCT_2";
        case DESCRIPTOR_PRODUCT_3:
            return "DESCRIPTOR_PRODUCT_3";
        case DESCRIPTOR_DETERMINANT:
            return "DESCRIPTOR_DETERMINANT";
        case DESCRIPTOR_CD:
            return "DESCRIPTOR_CD";
        case DESCRIPTOR_CUV:
            return "DESCRIPTOR_CUV";
        case DESCRIPTOR_CRK:
            return "DESCRIPTOR_CRK";
        case DESCRIPTOR_CK:
            return "DESCRIPTOR_CK";
        case DESCRIPTOR_CJ:
            return "DESCRIPTOR_CJ";
        case DESCRIPTOR_Z_OBSERVATIONS:
            return "DESCRIPTOR_Z_OBSERVATIONS";
        case DESCRIPTOR_Z_Actual:
            return "DESCRIPTOR_Z_Actual";
        case DESCRIPTOR_MSPE:
            return "DESCRIPTOR_MSPE";
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
        case DESCRIPTOR_MSPE_1:
            return "DESCRIPTOR_MSPE_1";
        case DESCRIPTOR_MSPE_2:
            return "DESCRIPTOR_MSPE_2";
        case DESCRIPTOR_Z_MISS:
            return "DESCRIPTOR_Z_MISS";
        case DESCRIPTOR_k_T :
            return "DESCRIPTOR_k_T";
        case DESCRIPTOR_k_A :
            return "DESCRIPTOR_k_A";
        case DESCRIPTOR_k_A_TMP :
            return "DESCRIPTOR_k_A_TMP";
        case DESCRIPTOR_k_T_TMP :
            return "DESCRIPTOR_k_T_TMP";
        case DESCRIPTOR_K_T :
            return "DESCRIPTOR_K_T";
        case DESCRIPTOR_K_T_TMP :
            return "DESCRIPTOR_K_T_TMP";
        case DESCRIPTOR_K_A :
            return "DESCRIPTOR_K_A";
        case DESCRIPTOR_EXPR_1 :
            return "DESCRIPTOR_EXPR_1";
        case DESCRIPTOR_EXPR_2 :
            return "DESCRIPTOR_EXPR_2";
        case DESCRIPTOR_EXPR_3 :
            return "DESCRIPTOR_EXPR_3";
        case DESCRIPTOR_EXPR_4 :
            return "DESCRIPTOR_EXPR_4";
        case DESCRIPTOR_MLOE :
            return "DESCRIPTOR_MLOE";
        case DESCRIPTOR_MMOM:
            return "DESCRIPTOR_MMOM";
        case DESCRIPTOR_ALPHA :
            return "DESCRIPTOR_ALPHA";
        case DESCRIPTOR_TRUTH_ALPHA :
            return "DESCRIPTOR_TRUTH_ALPHA";
        case DESCRIPTOR_TIMATED_ALPHA :
            return "DESCRIPTOR_TIMATED_ALPHA";
        case DESCRIPTOR_MLOE_MMOM :
            return "DESCRIPTOR_MLOE_MMOM";
        case DESCRIPTOR_A :
            return "DESCRIPTOR_A";
        case DESCRIPTOR_RESULTS :
            return "DESCRIPTOR_RESULTS";
        case DESCRIPTOR_C_TRACE :
            return "DESCRIPTOR_C_TRACE";
        case DESCRIPTOR_C_DIAG :
            return "DESCRIPTOR_C_DIAG";
        case DESCRIPTOR_SUM :
            return "DESCRIPTOR_SUM";
        case DESCRIPTOR_R :
            return "DESCRIPTOR_R";
        case DESCRIPTOR_R_COPY :
            return "DESCRIPTOR_R_COPY";
        case DESCRIPTOR_F_DATA:
            return "DESCRIPTOR_F_DATA";
        case DESCRIPTOR_ET1:
            return "DESCRIPTOR_ET1";
        case DESCRIPTOR_ET2:
            return "DESCRIPTOR_ET2";
        case DESCRIPTOR_EP:
            return "DESCRIPTOR_EP";
        case DESCRIPTOR_SLMN:
            return "DESCRIPTOR_SLMN";
        case DESCRIPTOR_IE:
            return "DESCRIPTOR_IE";
        case DESCRIPTOR_IO:
            return "DESCRIPTOR_IO";
        case DESCRIPTOR_P:
            return "DESCRIPTOR_P";
        case DESCRIPTOR_D:
            return "DESCRIPTOR_D";
        case DESCRIPTOR_FLMERA:
            return "DESCRIPTOR_FLMERA";
        case DESCRIPTOR_ZLM:
            return "DESCRIPTOR_ZLM";
        case DESCRIPTOR_SC:
            return "DESCRIPTOR_SC";
        case DESCRIPTOR_F_SPATIAL:
            return "DESCRIPTOR_F_SPATIAL";
        case DESCRIPTOR_FLM:
            return "DESCRIPTOR_FLM";
        case DESCRIPTOR_FLMT:
            return "DESCRIPTOR_FLMT";
        default:
            throw std::invalid_argument(
                    "The name of descriptor you provided is undefined, Please read the user manual to know the available descriptors");
    }
}

template<typename T>
bool DescriptorData<T>::GetIsDescriptorInitiated() {
    return this->mIsDescriptorInitiated;
}

template<typename T>
void DescriptorData<T>::SetIsDescriptorInitiated(bool aIsInitiated) {
    this->mIsDescriptorInitiated = aIsInitiated;
}