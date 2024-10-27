
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file DataTransformer.cpp
 * @brief Contains the implementation of the DataTransformer class.
 * @version 2.0.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2024-02-04
**/

#include <data-transformer/DataTransformer.hpp>
extern "C"{
    #include <runtime/parsec/jdf/JobDescriptionFormat.h>
}

using namespace exageostat::common;
using namespace exageostat::configurations;
using namespace exageostat::transformers;

template<typename T> void DataTransformer<T>::ForwardSphericalHarmonicsTransform(const int &aLSize, std::unique_ptr<ExaGeoStatData<T>> &aData){

    SYNC_TIME_START();
    parsec_tiled_matrix_t *pFDataDesc = (parsec_tiled_matrix_t *) aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_F_DATA).parsec_desc;
    parsec_tiled_matrix_t *pFLMDesc = (parsec_tiled_matrix_t *) aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_FLM).parsec_desc;
    parsec_tiled_matrix_t *pFLMTDesc = (parsec_tiled_matrix_t *) aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_FLMT).parsec_desc;
    parsec_tiled_matrix_t *pET1Desc = (parsec_tiled_matrix_t *) aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_ET1).parsec_desc;
    parsec_tiled_matrix_t *pET2Desc = (parsec_tiled_matrix_t *) aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_ET2).parsec_desc;
    parsec_tiled_matrix_t *pEPDesc = (parsec_tiled_matrix_t *) aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_EP).parsec_desc;
    parsec_tiled_matrix_t *pSLMNDesc = (parsec_tiled_matrix_t *) aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_SLMN).parsec_desc;
    parsec_tiled_matrix_t *pIEDesc = (parsec_tiled_matrix_t *) aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_IE).parsec_desc;
    parsec_tiled_matrix_t *pIODesc = (parsec_tiled_matrix_t *) aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_IO).parsec_desc;
    parsec_tiled_matrix_t *pPDesc = (parsec_tiled_matrix_t *) aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_P).parsec_desc;
    parsec_tiled_matrix_t *pDDesc = (parsec_tiled_matrix_t *) aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_D).parsec_desc;

    int aFDataM = pFDataDesc->mb;
    int aEPN = pEPDesc->nb;
    int aET1M = pET1Desc->mb;
    int aET2M = pET2Desc->mb;
    int aPN = pPDesc->nb;
    int aFlmM = pFLMDesc->mb;
    int aFlmN = pFLMDesc->nb;

    ForwardSHT((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), pFDataDesc, pFLMDesc,
               pFLMTDesc, pET1Desc, pET2Desc, pEPDesc, pSLMNDesc, pIEDesc, pIODesc, pPDesc,
               pDDesc, aFDataM, aEPN, aET1M, aET2M, aPN, aFlmM, aFlmN, aLSize);

}

template<typename T> void DataTransformer<T>::ForwardReshape(Configurations &aConfigurations, std::unique_ptr<ExaGeoStatData<T>> &aData){

    SYNC_TIME_START();
    int rank = ExaGeoStatHardware::GetParsecMPIRank();
    int verbose = configurations::Configurations::GetVerbosity() == DETAILED_MODE? 1: 0;
    int L = aConfigurations.GetDenseTileSize();
    int aT = aConfigurations.GetTimeSlot();

    parsec_tiled_matrix_t *pFDataDesc = (parsec_tiled_matrix_t *) aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_F_DATA).parsec_desc;
    parsec_tiled_matrix_t *pFLMDesc = (parsec_tiled_matrix_t *) aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_FLM).parsec_desc;
    parsec_tiled_matrix_t *pFLMTDesc = (parsec_tiled_matrix_t *) aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_FLMT).parsec_desc;
    parsec_tiled_matrix_t *pET1Desc = (parsec_tiled_matrix_t *) aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_ET1).parsec_desc;
    parsec_tiled_matrix_t *pET2Desc = (parsec_tiled_matrix_t *) aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_ET2).parsec_desc;
    parsec_tiled_matrix_t *pEPDesc = (parsec_tiled_matrix_t *) aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_EP).parsec_desc;
    parsec_tiled_matrix_t *pSLMNDesc = (parsec_tiled_matrix_t *) aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_SLMN).parsec_desc;
    parsec_tiled_matrix_t *pIEDesc = (parsec_tiled_matrix_t *) aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_IE).parsec_desc;
    parsec_tiled_matrix_t *pIODesc = (parsec_tiled_matrix_t *) aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_IO).parsec_desc;
    parsec_tiled_matrix_t *pPDesc = (parsec_tiled_matrix_t *) aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_P).parsec_desc;
    parsec_tiled_matrix_t *pDDesc = (parsec_tiled_matrix_t *) aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_D).parsec_desc;
    parsec_tiled_matrix_t *pADesc = (parsec_tiled_matrix_t *) aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_A).parsec_desc;

    int aFDataM = pFDataDesc->mb;
    int aEPN = pEPDesc->nb;
    int aET1M = pET1Desc->mb;
    int aET2M = pET2Desc->mb;
    int aPN = pPDesc->nb;
    int aFlmTM = pFLMTDesc->mb;
    int aFlmTN = pFLMTDesc->nb;
    int nodes = aConfigurations.GetCoresNumber();
    int aFlmTNB = (aFlmTM / nodes % L) ? aFlmTM / nodes / L + 1 : aFlmTM / nodes / L;

    double aNormGlobal = 0;
    double aNT = 0;
    double aUpperLower = EXAGEOSTAT_LOWER;
    ForwardSHTReshape((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), rank, verbose, pFDataDesc, pFLMDesc,
                      pFLMTDesc, pET1Desc, pET2Desc, pEPDesc, pSLMNDesc, pIEDesc, pIODesc, pPDesc,
                      pDDesc, pADesc, aFDataM, aEPN, aET1M, aET2M, aPN, aFlmTNB, aT, L, aNormGlobal, aNT, aUpperLower);
}

template<typename T> void DataTransformer<T>::InverseSphericalHarmonicsTransform(const int &aLSize, std::unique_ptr<ExaGeoStatData<T>> &aData){

    parsec_tiled_matrix_t *pFSpatialDesc = (parsec_tiled_matrix_t *) aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_F_SPATIAL).parsec_desc;
    parsec_tiled_matrix_t *pFLMDesc = (parsec_tiled_matrix_t *) aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_FLM).parsec_desc;
    parsec_tiled_matrix_t *pZLMDesc = (parsec_tiled_matrix_t *) aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_ZLM).parsec_desc;
    parsec_tiled_matrix_t *pSCDesc = (parsec_tiled_matrix_t *) aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_SC).parsec_desc;

    InverseSHT((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), pFSpatialDesc, pFLMDesc, pZLMDesc, pSCDesc, aLSize);

}
