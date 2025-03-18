
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file DataAnalyzer.cpp
 * @brief Contains the implementation of the DataAnalyzer class.
 * @version 2.0.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @author Qinglei Cao
 * @date 2024-02-04
**/

#include <data-analyzer/DataAnalyzer.hpp>

extern "C" {
#include <runtime/parsec/jdf/JobDescriptionFormat.h>
}

using namespace exageostat::analyzer;
using namespace exageostat::common;

template<typename T>
void
DataAnalyzer<T>::PreAnalyzeMatrix(std::unique_ptr<ExaGeoStatData<T>> &aData){

    auto *pParams = ExaGeoStatHardware::GetHicmaParams();
    auto *pHicma_data = ExaGeoStatHardware::GetHicmaData();
    auto *pAnalysis = ExaGeoStatHardware::GetAnalysis();
    auto *pStarsh_kernel = ExaGeoStatHardware::GetParamsKernel();
    auto *pContext = (parsec_context_t *) ExaGeoStatHardware::GetParsecContext();

    hicma_parsec_matrix_pre_analysis(pContext, pHicma_data, pParams, pStarsh_kernel, pAnalysis);
}

template<typename T>
void
DataAnalyzer<T>::PostAnalyzeMatrix(std::unique_ptr<ExaGeoStatData<T>> &aData){

    auto *pParams = ExaGeoStatHardware::GetHicmaParams();
    auto *pHicma_data = ExaGeoStatHardware::GetHicmaData();
    auto *pAnalysis = ExaGeoStatHardware::GetAnalysis();
    auto *pStarsh_kernel = ExaGeoStatHardware::GetParamsKernel();
    auto *pContext = (parsec_context_t *) ExaGeoStatHardware::GetParsecContext();

    hicma_parsec_matrix_post_analysis(pContext, pHicma_data, pParams, pStarsh_kernel, pAnalysis);
}

template<typename T>
double
DataAnalyzer<T>::CompareMatDifference(std::unique_ptr<ExaGeoStatData<T>> &aData){

    // Get parsec descriptors
    auto flm_desc = aData->GetDescriptorData()->GetDescriptor(DescriptorType::PARSEC_DESCRIPTOR,DescriptorName::DESCRIPTOR_FLM).parsec_desc;
    auto flm_era_desc = aData->GetDescriptorData()->GetDescriptor(DescriptorType::PARSEC_DESCRIPTOR,DescriptorName::DESCRIPTOR_FLMERA).parsec_desc;

    // Call jdf generated function
    DifferenceDouble((parsec_context_t *)ExaGeoStatHardware::GetParsecContext(), flm_desc, flm_era_desc);
}
