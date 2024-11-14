
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ParsecRuntimeSolver.cpp
 * @brief This file contains the implementation of ParsecRuntimeSolver class.
 * @details ParsecRuntimeSolver is a concrete implementation of the RuntimeSolversMethods class.
 * @version 2.0.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @author Qinglei Cao
 * @date 2024-11-04
**/

#include <runtime-solver/concrete/ParsecRuntimeSolver.hpp>
#include <data-analyzer/DataAnalyzer.hpp>
#include <data-transformer/DataTransformer.hpp>

extern "C"{
#include <runtime/parsec/jdf/JobDescriptionFormat.h>
}

using namespace exageostat::common;
using namespace exageostat::configurations;
using namespace exageostat::analyzer;
using namespace exageostat::runtimesolver;

template<typename T>
void ParsecRuntimeSolver<T>::ExaGeoStatSYRK(std::unique_ptr<ExaGeoStatData<T>> &aData){
    auto* pContext = (parsec_context_t *) ExaGeoStatHardware::GetParsecContext();
    auto* pDesc_A = (parsec_tiled_matrix_t *) aData->GetDescriptorData()->GetDescriptor(DescriptorType::PARSEC_DESCRIPTOR, DescriptorName::DESCRIPTOR_A).parsec_desc;

    SYNC_TIME_START();
    dplasma_dsyrk(pContext, dplasmaLower, dplasmaNoTrans, 1.0, pDesc_A, 0.0, (parsec_tiled_matrix_t *) &ExaGeoStatHardware::GetHicmaData()->dcA);
    SYNC_TIME_PRINT(ExaGeoStatHardware::GetParsecMPIRank(),("SYRK\n"));
}

template<typename T>
void ParsecRuntimeSolver<T>::ExaGeoStatTLRCholesky(std::unique_ptr<ExaGeoStatData<T>> &aData){

    auto *pParams = ExaGeoStatHardware::GetHicmaParams();
    auto *pHicma_data = ExaGeoStatHardware::GetHicmaData();
    auto *pAnalysis = ExaGeoStatHardware::GetAnalysis();
    auto *pContext = (parsec_context_t *) ExaGeoStatHardware::GetParsecContext();
    for( int i= 0; i < pParams->nruns; i++ ) {
        hicma_parsec_potrf(pContext, pHicma_data, pParams, pAnalysis);
    }
}

template<typename T>
double ParsecRuntimeSolver<T>::ExaGeoStatNorm(Configurations &aConfigurations, std::unique_ptr<ExaGeoStatData<T>> &aData){

    int L = aConfigurations.GetDenseTileSize();
    int N = aConfigurations.GetProblemSize();
    double aNT = (N % L == 0) ? (N/L) : (N/L + 1);
    int aUpperLower = EXAGEOSTAT_LOWER;
    auto* pContext = (parsec_context_t *) ExaGeoStatHardware::GetParsecContext();
    auto* pDescA = (parsec_tiled_matrix_t *) aData->GetDescriptorData()->GetDescriptor(DescriptorType::PARSEC_DESCRIPTOR, DescriptorName::DESCRIPTOR_A).parsec_desc;

    SYNC_TIME_START();
    GetMatrixNorm(pContext, &ExaGeoStatHardware::GetHicmaParams()->norm_global, (parsec_tiled_matrix_t *) &ExaGeoStatHardware::GetHicmaData()->dcA, aNT, aUpperLower, 1);
    SYNC_TIME_PRINT(ExaGeoStatHardware::GetParsecMPIRank(), ("Matrix norm: norm_global= %le\n", ExaGeoStatHardware::GetHicmaParams()->norm_global));
    return ExaGeoStatHardware::GetHicmaParams()->norm_global;
}

template<typename T>
double ParsecRuntimeSolver<T>::CalculateMSE(Configurations &aConfigurations, std::unique_ptr<ExaGeoStatData<T>> &aData) {

    auto* pContext = (parsec_context_t * )ExaGeoStatHardware::GetParsecContext();
    auto* pDesc_f_data = aData->GetDescriptorData()->GetDescriptor(DescriptorType::PARSEC_DESCRIPTOR,
                                                                  DescriptorName::DESCRIPTOR_F_DATA).parsec_desc;
    auto* pDesc_f_spatial = aData->GetDescriptorData()->GetDescriptor(DescriptorType::PARSEC_DESCRIPTOR,
                                                                     DescriptorName::DESCRIPTOR_F_SPATIAL).parsec_desc;

    SYNC_TIME_START();
    auto mse_result = MeanSquaredError(pContext, pDesc_f_data, pDesc_f_spatial, aConfigurations.GetDenseTileSize());
    SYNC_TIME_PRINT(ExaGeoStatHardware::GetParsecMPIRank(),("mse\n"));
    return mse_result;
}

template<typename T>
T ParsecRuntimeSolver<T>::ModelingOperations(std::unique_ptr<ExaGeoStatData<T>> &aData, Configurations &aConfigurations,
                                              T *apMeasurementsMatrix, const kernels::Kernel<T> &aKernel) {

     // SYRK
    ExaGeoStatSYRK(aData);
    // Calculate norm
    ExaGeoStatNorm(aConfigurations, aData);
    // Analyze matrix before Cholesky
    DataAnalyzer<T>::PreAnalyzeMatrix(aData);
    // HiCMA Cholesky
    ExaGeoStatTLRCholesky(aData);
    // Analyze matrix after Cholesky
    DataAnalyzer<T>::PostAnalyzeMatrix(aData);
    // Diff to matlab result
    DataAnalyzer<T>::CompareMatDifference(aData);

    if(aConfigurations.GetEnableInverse()){
        transformers::DataTransformer<T>::InverseSphericalHarmonicsTransform(aConfigurations.GetDenseTileSize(), aData);
        // TODO: results in a seg fault in C
        CalculateMSE(aConfigurations, aData);
    }
}
