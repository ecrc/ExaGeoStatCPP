
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ClimateEmulator.cpp
 * @brief High-Level Wrapper class containing the static API for ClimateEmulator operations.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @date 2025-01-19
**/

#include <api/ClimateEmulator.hpp>
#include <prediction/Prediction.hpp>
#include <data-analyzer/DataAnalyzer.hpp>
#include <data-transformer/DataTransformer.hpp>
#include <data-loader/concrete/ParsecLoader.hpp>
#include <runtime-solver/concrete/ParsecRuntimeSolver.hpp>

extern "C" {
#include <runtime/parsec/jdf/JobDescriptionFormat.h>
}

#ifdef USE_CUDA
#include <runtime/parsec/GPUHelperFunctions.h>
#endif

using namespace std;

using namespace exageostat::api;
using namespace exageostat::common;
using namespace exageostat::analyzer;
using namespace exageostat::dataunits;
using namespace exageostat::dataLoader;
using namespace exageostat::transformers;
using namespace exageostat::runtimesolver;
using namespace exageostat::configurations;

template<typename T>
void ClimateEmulator<T>::ClimateEmulatorLoadData(Configurations &aConfigurations,
                                                 unique_ptr <ExaGeoStatData<T>> &aData) {

    aConfigurations.PrintSummary();
    LOGGER("** Climate Emulator data generation/loading **")

    // Create data object
    if (!aData) {
        aData = std::make_unique < ExaGeoStatData <
                T >> (aConfigurations.GetProblemSize() / 1, aConfigurations.GetDimension());
    }

    SYNC_TIME_START();

    // Initiate Descriptors
    int L = aConfigurations.GetDenseTileSize();
    int MB;
    int NB;
    int t = aConfigurations.GetTimeSlot();
    int P = aConfigurations.GetPGrid();
    int nodes = aConfigurations.GetCoresNumber();
    int rank = ExaGeoStatHardware::GetParsecMPIRank();
    int verbose = configurations::Configurations::GetVerbosity() == DETAILED_MODE ? 1 : 0;
    int gpus = aConfigurations.GetGPUsNumbers();
    int tile_size = aConfigurations.GetDenseTileSize() * 10;
    string files_directory_path = aConfigurations.GetDataPath();
    int path_length = files_directory_path.length();
    char filename[path_length + 50];
    char directory_path[path_length];
    sprintf(directory_path, "%s", files_directory_path.c_str());

    MB = L + 1;
    NB = L * 2;
    VERBOSE_PRINT(rank, verbose, ("Reading f_data\n"));
    aData->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_F_DATA);
    parsec_matrix_block_cyclic_t *pF_data_desc = aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR,
                                                                                           DESCRIPTOR_F_DATA).parsec_desc;
    ReadCSVToComplexTimeSlot((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), pF_data_desc, MB, NB, nodes,
                             t, directory_path, rank, verbose, gpus);

    MB = 2 * L - 1;
    NB = L + 1;
    VERBOSE_PRINT(rank, verbose, ("Reading Et1\n"));
    aData->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_ET1);
    parsec_matrix_block_cyclic_t *pEt1_data_desc = aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR,
                                                                                             DESCRIPTOR_ET1).parsec_desc;
    sprintf(filename, "%s/%d%s", files_directory_path.c_str(), tile_size, "_Et1.csv");
    ReadCSVComplex((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), pEt1_data_desc, MB, NB, nodes, t,
                   filename, rank, verbose, gpus);

    MB = 2 * L - 1;
    NB = L - 1;
    VERBOSE_PRINT(rank, verbose, ("Reading Et2\n"));
    aData->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_ET2);
    parsec_matrix_block_cyclic_t *pEt2_data_desc = aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR,
                                                                                             DESCRIPTOR_ET2).parsec_desc;
    sprintf(filename, "%s/%d%s", files_directory_path.c_str(), tile_size, "_Et2.csv");
    ReadCSVComplex((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), pEt2_data_desc, MB, NB, nodes, t,
                   filename, rank, verbose, gpus);

    MB = 2 * L;
    NB = 2 * L - 1;
    VERBOSE_PRINT(rank, verbose, ("Reading Ep\n"));
    aData->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_EP);
    parsec_matrix_block_cyclic_t *pEp_data_desc = aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR,
                                                                                            DESCRIPTOR_EP).parsec_desc;
    sprintf(filename, "%s/%d%s", files_directory_path.c_str(), tile_size, "_Ep.csv");
    ReadCSVComplex((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), pEp_data_desc, MB, NB, nodes, t,
                   filename, rank, verbose, gpus);

    MB = (L * L + L) / 2;
    NB = L;
    VERBOSE_PRINT(rank, verbose, ("Reading Slmn\n"));
    aData->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_SLMN);
    parsec_matrix_block_cyclic_t *pSlum_data_desc = aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR,
                                                                                              DESCRIPTOR_SLMN).parsec_desc;
    sprintf(filename, "%s/%d%s", files_directory_path.c_str(), tile_size, "_Slmn.csv");
    ReadCSVComplex((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), pSlum_data_desc, MB, NB, nodes, t,
                   filename, rank, verbose, gpus);

    MB = L;
    NB = 2 * L - 1;
    VERBOSE_PRINT(rank, verbose, ("Reading Ie\n"));
    aData->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_IE);
    parsec_matrix_block_cyclic_t *pIe_data_desc = aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR,
                                                                                            DESCRIPTOR_IE).parsec_desc;
    sprintf(filename, "%s/%d%s", files_directory_path.c_str(), tile_size, "_Ie.csv");
    ReadCSVToComplex((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), pIe_data_desc, MB, NB, nodes, t,
                     filename, rank, verbose, gpus);

    MB = L;
    NB = 2 * L - 1;
    VERBOSE_PRINT(rank, verbose, ("Reading Io\n"));
    aData->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_IO);
    parsec_matrix_block_cyclic_t *pIo_data_desc = aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR,
                                                                                            DESCRIPTOR_IO).parsec_desc;
    sprintf(filename, "%s/%d%s", files_directory_path.c_str(), tile_size, "_Io.csv");
    ReadCSVToComplex((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), pIo_data_desc, MB, NB, nodes, t,
                     filename, rank, verbose, gpus);

    MB = L - 1;
    NB = L + 1;
    VERBOSE_PRINT(rank, verbose, ("Reading P\n"));
    aData->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_P);
    parsec_matrix_block_cyclic_t *pP_data_desc = aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR,
                                                                                           DESCRIPTOR_P).parsec_desc;
    sprintf(filename, "%s/%d%s", files_directory_path.c_str(), tile_size, "_P.csv");
    ReadCSVToComplex((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), pP_data_desc, MB, NB, nodes, t,
                     filename, rank, verbose, gpus);

    MB = 2 * L - 1;
    NB = 2 * L - 1;
    VERBOSE_PRINT(rank, verbose, ("Reading D\n"));
    aData->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_D);
    parsec_matrix_block_cyclic_t *pD_data_desc = aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR,
                                                                                           DESCRIPTOR_D).parsec_desc;
    sprintf(filename, "%s/%d%s", files_directory_path.c_str(), tile_size, "_D.csv");
    ReadCSVToComplex((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), pD_data_desc, MB, NB, nodes, t,
                     filename, rank, verbose, gpus);

    MB = L;
    NB = L;
    VERBOSE_PRINT(rank, verbose, ("Reading flm\n"));
    aData->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_FLM);
    parsec_matrix_block_cyclic_t *pFlm_data_desc = aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR,
                                                                                             DESCRIPTOR_FLM).parsec_desc;
    ReadCSVTimeSlot((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), pFlm_data_desc, MB, NB, nodes, t,
                    filename, rank, verbose, gpus);

    VERBOSE_PRINT(rank, verbose, ("Reading flmERA\n"));
    aData->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_FLMERA);
    parsec_matrix_block_cyclic_t *pFlmera_data_desc = aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR,
                                                                                                DESCRIPTOR_FLMERA).parsec_desc;
    sprintf(filename, "%s/%d%s", files_directory_path.c_str(), tile_size, "_flmERA.csv");
    ReadCSVTimeSlot((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), pFlmera_data_desc, MB, NB, nodes, t,
                    filename, rank, verbose, gpus);

    // Backward
    if (aConfigurations.GetEnableInverse()) {

        MB = L + 1;
        NB = (L * L + L) / 2;
        VERBOSE_PRINT(rank, verbose, ("Reading Zlm\n"));
        aData->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_ZLM);
        parsec_matrix_block_cyclic_t *PZlm_data_desc = aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR,
                                                                                                 DESCRIPTOR_ZLM).parsec_desc;
        sprintf(filename, "%s/%d%s", files_directory_path.c_str(), tile_size, "_Zlm.csv");
        ReadCSV((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), PZlm_data_desc, MB, NB, nodes, t, filename,
                rank, verbose, gpus);

        MB = 2 * L - 1;
        NB = 2 * L;
        VERBOSE_PRINT(rank, verbose, ("Reading SC\n"));
        aData->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_SC);
        parsec_matrix_block_cyclic_t *pSc_data_desc = aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR,
                                                                                                DESCRIPTOR_SC).parsec_desc;
        sprintf(filename, "%s/%d%s", files_directory_path.c_str(), tile_size, "_SC.csv");
        ReadCSV((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), pSc_data_desc, MB, NB, nodes, t, filename,
                rank, verbose, gpus);

        MB = L + 1;
        NB = 2 * L;
        VERBOSE_PRINT(rank, verbose, ("f_spatial\n"));
        aData->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_F_SPATIAL);
        parsec_matrix_block_cyclic_t *pF_spatial_data_desc = aData->GetDescriptorData()->GetDescriptor(
                PARSEC_DESCRIPTOR, DESCRIPTOR_F_SPATIAL).parsec_desc;
        ReadCSV((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), pF_spatial_data_desc, MB, NB, nodes, t,
                filename, rank, verbose, gpus);
    }

    // Init and allocate memory for desc_flmT
    MB = L * L;
    NB = t;
    aData->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_FLMT);
    parsec_matrix_block_cyclic_t *pFlmt_data_desc = aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR,
                                                                                              DESCRIPTOR_FLMT).parsec_desc;
    parsec_matrix_block_cyclic_init(pFlmt_data_desc, PARSEC_MATRIX_DOUBLE, PARSEC_MATRIX_TILE, rank,
                                    L * ((MB / nodes % L) ? MB / nodes / L + 1 : MB / nodes / L),
                                    NB, MB, NB, 0, 0, MB, NB, nodes, 1, 1, 1, 0, 0);

    pFlmt_data_desc->mat = parsec_data_allocate((size_t) pFlmt_data_desc->super.nb_local_tiles *
                                                (size_t) pFlmt_data_desc->super.bsiz *
                                                (size_t) parsec_datadist_getsizeoftype(pFlmt_data_desc->super.mtype));
    parsec_data_collection_set_key((parsec_data_collection_t * ) & desc_flmT, "desc_flmT");

    // Init and allocate memory for pA_data_desc
    MB = L * L;
    NB = t;
    aData->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_A);
    parsec_matrix_block_cyclic_t *pA_data_desc = aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR,
                                                                                           DESCRIPTOR_A).parsec_desc;
    parsec_matrix_block_cyclic_init(pA_data_desc, PARSEC_MATRIX_DOUBLE, PARSEC_MATRIX_TILE, rank, L, L, MB, NB, 0, 0,
                                    pFlmt_data_desc->super.mb, pFlmt_data_desc->super.nb, P, nodes / P, 1, 1, 0, 0);
    pA_data_desc->mat = parsec_data_allocate((size_t) pA_data_desc->super.nb_local_tiles *
                                             (size_t) pA_data_desc->super.bsiz *
                                             (size_t) parsec_datadist_getsizeoftype(pA_data_desc->super.mtype));
    parsec_data_collection_set_key((parsec_data_collection_t * ) & pA_data_desc, "desc_A");

    if (aConfigurations.GetEnableInverse()) {
        int ts_test_M = 2000;
        int ts_test_N = 1;
        // Allocate memory
        double *pFileContent = (double *) malloc(ts_test_M * ts_test_N * sizeof(double));
        sprintf(filename, "%s/%s", files_directory_path.c_str(), "ts_test.csv");
        parsec::ParsecLoader<T>::ReadCSVFileHelper(filename, pFileContent, ts_test_M,
                                                   ts_test_N);
    }

    SYNC_TIME_PRINT(ExaGeoStatHardware::GetParsecMPIRank(), ("Load Data\n"));
    void *work_space = nullptr;
#ifdef USE_CUDA

    if(aConfigurations.GetGPUsNumbers() > 0) {
        work_space = InitializeWorkSpace((parsec_context_t *)ExaGeoStatHardware::GetParsecContext(), pF_data_desc->super.mb,
                                              pEp_data_desc->super.nb, pEt1_data_desc->super.mb, pEt2_data_desc->super.mb,
                                              pP_data_desc->super.nb);
    }
#endif

    // Forward SHT
    DataTransformer<T>::ForwardSphericalHarmonicsTransform(aConfigurations, aData, work_space);
    // Forward SHT Reshape
    DataTransformer<T>::ForwardReshape(aConfigurations, aData);
    // Generate matrix
    parsec::ParsecLoader<T>::CompressMatrixHelper(aConfigurations, aData);

    LOGGER("\t*Data generation/loading finished*")
}

template<typename T>
void
ClimateEmulator<T>::ClimateEmulatorModeling(Configurations &aConfigurations, unique_ptr <ExaGeoStatData<T>> &aData,
                                            T *apMeasurementsMatrix) {

    aConfigurations.PrintSummary();
    LOGGER("** Climate Emulator data Modeling **")

    // Initialize all theta: starting, estimated, lower and upper bounds.
    aConfigurations.InitializeAllTheta();

    ParsecRuntimeSolver<T>::ExaGeoStatSYRK(aData);
    // Calculate norm
    ParsecRuntimeSolver<T>::ExaGeoStatNorm(aConfigurations, aData);
    // Analyze matrix before Cholesky
    DataAnalyzer<T>::PreAnalyzeMatrix(aData);
    // HiCMA Cholesky
    ParsecRuntimeSolver<T>::ExaGeoStatTLRCholesky(aData);
    // Analyze matrix after Cholesky
    DataAnalyzer<T>::PostAnalyzeMatrix(aData);
    // Diff to matlab result
    DataAnalyzer<T>::CompareMatDifference(aData);

    if (aConfigurations.GetEnableInverse()) {
        transformers::DataTransformer<T>::InverseSphericalHarmonicsTransform(aConfigurations.GetDenseTileSize(), aData);
        // TODO: results in a seg fault in C
        ParsecRuntimeSolver<T>::CalculateMSE(aConfigurations, aData);
    }
    LOGGER("\t*Data modeling finished*")
}


