
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file DataTransformer.cpp
 * @brief Contains the implementation of the DataTransformer class.
 * @version 2.0.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @author Qinglei Cao
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

    int f_desc_M = pFDataDesc->mb;
    int ep_desc_N = pEPDesc->nb;
    int et1_desc_M = pET1Desc->mb;
    int et2_desc_M = pET2Desc->mb;
    int p_desc_N = pPDesc->nb;
    int flm_desc_M = pFLMDesc->mb;
    int flm_desc_N = pFLMDesc->nb;

    ForwardSHT((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), pFDataDesc, pFLMDesc,
               pFLMTDesc, pET1Desc, pET2Desc, pEPDesc, pSLMNDesc, pIEDesc, pIODesc, pPDesc,
               pDDesc, f_desc_M, ep_desc_N, et1_desc_M, et2_desc_M, p_desc_N, flm_desc_M, flm_desc_N, aLSize);

    int flops_forward = 2.0*(aLSize+1)*(2*aLSize-1)*(2*aLSize) // Gmtheta_r = f_data*Ep
        + 2.0*(2*aLSize-1)*(2*aLSize-1)*(aLSize+1) // Fmnm = Et1*Gmtheta_r
        + 2.0*(2*aLSize-1)*(aLSize-1)*(aLSize+1) // tmp1 = Et2*P
        + 2.0*(2*aLSize-1)*(2*aLSize-1)*(aLSize+1)  // tmp2 = tmp1 * Gmtheta_r
        + 2.0*(2*aLSize-1)*(2*aLSize-1)*(2*aLSize-1)  // Fmnm += tmp2 * D
        + 2.0*aLSize*aLSize/2*(aLSize*(2*aLSize-1)+aLSize);   // flmn_matrix(ell+1,m+1) = Slmn(climate_emulator_getSingleIndex(ell, m),:)*Ie*Fmnm(:,L+m)

    SYNC_TIME_PRINT(ExaGeoStatHardware::GetParsecMPIRank(), ("Forward SHT: %.2lf Gflop/s\n", flops_forward/sync_time_elapsed/1.0e9));

}

template<typename T> void DataTransformer<T>::ForwardReshape(Configurations &aConfigurations, std::unique_ptr<ExaGeoStatData<T>> &aData){

    SYNC_TIME_START();
    int rank = ExaGeoStatHardware::GetParsecMPIRank();
    int verbose = configurations::Configurations::GetVerbosity() == DETAILED_MODE? 1: 0;
    int L = aConfigurations.GetDenseTileSize();
    int t = aConfigurations.GetTimeSlot();

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

    int f_desc_M = pFDataDesc->mb;
    int ep_desc_N = pEPDesc->nb;
    int et1_desc_M = pET1Desc->mb;
    int et2_desc_M = pET2Desc->mb;
    int p_desc_N = pPDesc->nb;
    int flmt_desc_M = pFLMTDesc->mb;
    int nodes = aConfigurations.GetCoresNumber();
    int flmt_desc_nb = (flmt_desc_M / nodes % L) ? flmt_desc_M / nodes / L + 1 : flmt_desc_M / nodes / L;

    int N = aConfigurations.GetProblemSize();
    int NT = (N % L == 0) ? (N/L) : (N/L + 1);
    double upper_lower = EXAGEOSTAT_LOWER;
    ForwardSHTReshape((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), rank, verbose, pFDataDesc, pFLMDesc,
                      pFLMTDesc, pET1Desc, pET2Desc, pEPDesc, pSLMNDesc, pIEDesc, pIODesc, pPDesc,
                      pDDesc, pADesc, f_desc_M, ep_desc_N, et1_desc_M, et2_desc_M, p_desc_N, flmt_desc_nb, t, L, &ExaGeoStatHardware::GetHicmaParams()->norm_global, NT, upper_lower);
    SYNC_TIME_PRINT(ExaGeoStatHardware::GetParsecMPIRank(), ("Forward SHT Reshape\n"));
    
}

template<typename T> void DataTransformer<T>::InverseSphericalHarmonicsTransform(const int &aLSize, std::unique_ptr<ExaGeoStatData<T>> &aData){

    SYNC_TIME_START();
    parsec_tiled_matrix_t *pFSpatialDesc = (parsec_tiled_matrix_t *) aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_F_SPATIAL).parsec_desc;
    parsec_tiled_matrix_t *pFLMDesc = (parsec_tiled_matrix_t *) aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_FLM).parsec_desc;
    parsec_tiled_matrix_t *pZLMDesc = (parsec_tiled_matrix_t *) aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_ZLM).parsec_desc;
    parsec_tiled_matrix_t *pSCDesc = (parsec_tiled_matrix_t *) aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_SC).parsec_desc;

    InverseSHT((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), pFSpatialDesc, pFLMDesc, pZLMDesc, pSCDesc, aLSize);
    SYNC_TIME_PRINT(ExaGeoStatHardware::GetParsecMPIRank(), ("Inverse SHT\n"));
}
