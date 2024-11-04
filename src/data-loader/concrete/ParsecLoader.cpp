
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file SyntheticGenerator.cpp
 * @brief Implementation of the SyntheticGenerator class
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @author Qinglei Cao
 * @date 2023-02-14
**/

#include <data-transformer/DataTransformer.hpp>
#include <data-loader/concrete/ParsecLoader.hpp>
extern "C"{
#include <runtime/parsec/jdf/JobDescriptionFormat.h>
}

using namespace std;

using namespace exageostat::dataLoader::parsec;
using namespace exageostat::common;
using namespace exageostat::configurations;
using namespace exageostat::dataunits;
using namespace exageostat::transformers;

template<typename T>
ParsecLoader<T> *ParsecLoader<T>::GetInstance() {

    if (mpInstance == nullptr) {
        mpInstance = new ParsecLoader<T>();
    }
    return mpInstance;
}

template<typename T>
std::unique_ptr<ExaGeoStatData<T>>
ParsecLoader<T>::LoadData(configurations::Configurations &aConfigurations, exageostat::kernels::Kernel<T> &aKernel) {

    SYNC_TIME_START();
    //create data object
    auto data = std::make_unique<ExaGeoStatData<T>>(aConfigurations.GetProblemSize() / 1,
                                                    aConfigurations.GetDimension());

    // Initiate Descriptors
    int L = aConfigurations.GetDenseTileSize();
    int MB;
    int NB;
    int t = aConfigurations.GetTimeSlot();
    int P = aConfigurations.GetPGrid();
    int nodes = aConfigurations.GetCoresNumber();
    int rank = ExaGeoStatHardware::GetParsecMPIRank();
    int verbose = configurations::Configurations::GetVerbosity() == DETAILED_MODE? 1: 0;
    int gpus = aConfigurations.GetGPUsNumbers();
    string files_directory_path = aConfigurations.GetDataPath();
    int path_length = files_directory_path.length();
    char filename[path_length + 50];
    char directory_path[path_length];
    sprintf(directory_path, "%s", files_directory_path.c_str());

    MB = L + 1;
    NB = L * 2;
	VERBOSE_PRINT(rank, verbose, ("Reading f_data\n"));
    data->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_F_DATA);
    parsec_matrix_block_cyclic_t *pF_data_desc = data->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_F_DATA).parsec_desc;
    ReadCSVToComplexTimeSlot((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), pF_data_desc, MB, NB, nodes, t, directory_path, rank, verbose, gpus);

    MB = 2*L-1;
    NB = L+1;
	VERBOSE_PRINT(rank, verbose, ("Reading Et1\n"));
    data->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_ET1);
    parsec_matrix_block_cyclic_t *pEt1_data_desc = data->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_ET1).parsec_desc;
    sprintf(filename, "%s/%s", files_directory_path.c_str(), "720_Et1.csv");
    ReadCSVComplex((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), pEt1_data_desc, MB, NB, nodes, t, filename, rank, verbose, gpus);

    MB = 2*L-1;
    NB = L-1;
	VERBOSE_PRINT(rank, verbose, ("Reading Et2\n"));
    data->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_ET2);
    parsec_matrix_block_cyclic_t *pEt2_data_desc = data->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_ET2).parsec_desc;
    sprintf(filename, "%s/%s", files_directory_path.c_str(), "720_Et2.csv");
    ReadCSVComplex((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), pEt2_data_desc, MB, NB, nodes, t, filename, rank, verbose, gpus);

    MB = 2*L;
    NB = 2*L-1;
	VERBOSE_PRINT(rank, verbose, ("Reading Ep\n"));
    data->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_EP);
    parsec_matrix_block_cyclic_t *pEp_data_desc = data->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_EP).parsec_desc;
    sprintf(filename, "%s/%s", files_directory_path.c_str(), "720_Ep.csv");
    ReadCSVComplex((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), pEp_data_desc, MB, NB, nodes, t, filename, rank, verbose, gpus);

    MB = (L*L+L)/2;
    NB = L;
	VERBOSE_PRINT(rank, verbose, ("Reading Slmn\n"));
    data->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_SLMN);
    parsec_matrix_block_cyclic_t *pSlum_data_desc = data->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_SLMN).parsec_desc;
    sprintf(filename, "%s/%s", files_directory_path.c_str(), "720_Slmn.csv");
    ReadCSVComplex((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), pSlum_data_desc, MB, NB, nodes, t, filename, rank, verbose, gpus);

    MB = L;
    NB = 2*L-1;
	VERBOSE_PRINT(rank, verbose, ("Reading Ie\n"));
    data->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_IE);
    parsec_matrix_block_cyclic_t *pIe_data_desc = data->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_IE).parsec_desc;
    sprintf(filename, "%s/%s", files_directory_path.c_str(), "720_Ie.csv");
    ReadCSVToComplex((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), pIe_data_desc, MB, NB, nodes, t, filename, rank, verbose, gpus);

    MB = L;
    NB = 2*L-1;
	VERBOSE_PRINT(rank, verbose, ("Reading Io\n"));
    data->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_IO);
    parsec_matrix_block_cyclic_t *pIo_data_desc = data->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_IO).parsec_desc;
    sprintf(filename, "%s/%s", files_directory_path.c_str(), "720_Io.csv");
    ReadCSVToComplex((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), pIo_data_desc, MB, NB, nodes, t, filename, rank, verbose, gpus);

    MB = L-1;
    NB = L+1;
	VERBOSE_PRINT(rank, verbose, ("Reading P\n"));
    data->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_P);
    parsec_matrix_block_cyclic_t *pP_data_desc = data->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_P).parsec_desc;
    sprintf(filename, "%s/%s", files_directory_path.c_str(), "720_P.csv");
    ReadCSVToComplex((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), pP_data_desc, MB, NB, nodes, t, filename, rank, verbose, gpus);

    MB = 2*L-1;
    NB = 2*L-1;
	VERBOSE_PRINT(rank, verbose, ("Reading D\n"));
    data->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_D);
    parsec_matrix_block_cyclic_t *pD_data_desc = data->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_D).parsec_desc;
    sprintf(filename, "%s/%s", files_directory_path.c_str(), "720_D.csv");
    ReadCSVToComplex((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), pD_data_desc, MB, NB, nodes, t, filename, rank, verbose, gpus);

    MB = L;
    NB = L;
	VERBOSE_PRINT(rank, verbose, ("Reading flm\n"));
    data->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_FLM);
    parsec_matrix_block_cyclic_t *pFlm_data_desc = data->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_FLM).parsec_desc;
    ReadCSVTimeSlot((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), pFlm_data_desc, MB, NB, nodes, t, filename, rank, verbose, gpus);

	VERBOSE_PRINT(rank, verbose, ("Reading flmERA\n"));
    data->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_FLMERA);
    parsec_matrix_block_cyclic_t *pFlmera_data_desc = data->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_FLMERA).parsec_desc;
    sprintf(filename, "%s/%s", files_directory_path.c_str(), "720_flmERA.csv");
    ReadCSVTimeSlot((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), pFlmera_data_desc, MB, NB, nodes, t, filename, rank, verbose, gpus);

    // Backward
    if(aConfigurations.GetEnableInverse()){

        MB = L+1;
        NB = (L*L+L)/2;
        VERBOSE_PRINT(rank, verbose, ("Reading Zlm\n"));
        data->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_ZLM);
        parsec_matrix_block_cyclic_t *PZlm_data_desc = data->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_ZLM).parsec_desc;
        sprintf(filename, "%s/%s", files_directory_path.c_str(), "720_Zlm.csv");
        ReadCSV((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), PZlm_data_desc, MB, NB, nodes, t, filename, rank, verbose, gpus);

        MB = 2*L-1;
        NB = 2*L;
        VERBOSE_PRINT(rank, verbose, ("Reading SC\n"));
        data->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_SC);
        parsec_matrix_block_cyclic_t *pSc_data_desc = data->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_SC).parsec_desc;
        sprintf(filename, "%s/%s", files_directory_path.c_str(), "720_SC.csv");
        ReadCSV((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), pSc_data_desc, MB, NB, nodes, t, filename, rank, verbose, gpus);

        MB = L+1;
        NB = 2*L;
        VERBOSE_PRINT(rank, verbose, ("f_spatial\n"));
        data->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_F_SPATIAL);
        parsec_matrix_block_cyclic_t *pF_spatial_data_desc = data->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_F_SPATIAL).parsec_desc;
        ReadCSV((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), pF_spatial_data_desc, MB, NB, nodes, t, filename, rank, verbose, gpus);
    }

    // Init and allocate memory for desc_flmT
    MB = L * L;
    NB = t;
    data->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_FLMT);
    parsec_matrix_block_cyclic_t *pFlmt_data_desc = data->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_FLMT).parsec_desc;
    parsec_matrix_block_cyclic_init(pFlmt_data_desc, PARSEC_MATRIX_DOUBLE, PARSEC_MATRIX_TILE, rank, L*((MB/nodes%L) ? MB/nodes/L+1 : MB/nodes/L),
                                    NB, MB, NB, 0, 0, MB, NB, nodes, 1, 1, 1, 0, 0);

    pFlmt_data_desc->mat = parsec_data_allocate((size_t)pFlmt_data_desc->super.nb_local_tiles *
                                   (size_t)pFlmt_data_desc->super.bsiz *
                                   (size_t)parsec_datadist_getsizeoftype(pFlmt_data_desc->super.mtype));
    parsec_data_collection_set_key((parsec_data_collection_t*)&desc_flmT, "desc_flmT");

    // Init and allocate memory for pA_data_desc
    MB = L * L;
    NB = t;
    data->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_A);
    parsec_matrix_block_cyclic_t *pA_data_desc = data->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_A).parsec_desc;
    parsec_matrix_block_cyclic_init(pA_data_desc, PARSEC_MATRIX_DOUBLE, PARSEC_MATRIX_TILE, rank, L, L, MB, NB, 0, 0,
                                    pFlmt_data_desc->super.mb, pFlmt_data_desc->super.nb, P, nodes/P, 1, 1, 0, 0);
    pA_data_desc->mat = parsec_data_allocate((size_t)pA_data_desc->super.nb_local_tiles *
                                   (size_t)pA_data_desc->super.bsiz *
                                   (size_t)parsec_datadist_getsizeoftype(pA_data_desc->super.mtype));
    parsec_data_collection_set_key((parsec_data_collection_t*)&pA_data_desc, "desc_A");

    if(aConfigurations.GetEnableInverse()){
        int ts_test_M = 2000;
        int ts_test_N = 1;
        // Allocate memory
        double *pFileContent = (double *)malloc(ts_test_M * ts_test_N * sizeof(double));
        sprintf(filename, "%s/%s", files_directory_path.c_str(), "ts_test.csv");
        ReadCSVFileHelper(filename, pFileContent, ts_test_M, ts_test_N);
    }

    SYNC_TIME_PRINT(ExaGeoStatHardware::GetParsecMPIRank(), ("Load Data\n"));
    // Forward SHT
    DataTransformer<T>::ForwardSphericalHarmonicsTransform(L, data);
    // Forward SHT Reshape
    DataTransformer<T>::ForwardReshape(aConfigurations, data);
    // Generate matrix
    CompressMatrixHelper(aConfigurations, data);

    return data;
}

template<typename T>
int ParsecLoader<T>::ReadCSVFileHelper(const char* apFilename, double *apFileContent, int aM, int aN) {

    FILE *pFile = fopen(apFilename, "r");
    if (!pFile) {
        printf("File opening failed: %s", apFilename);
        return -1;
    }

    int status = 0;
    for (int i = 0; i < aM; i++) {
        for (int j = 0; j < aN; j++) {
            status = fscanf(pFile, "%lf,", &apFileContent[i * aN + j]);
            if (status != 1) {
                fprintf(stderr, "Error reading file at row %d, column %d\n", i, j);
                fclose(pFile);
                return 1;
            }
        }
    }
    fclose(pFile);

    return 0;
}

template<typename T>
void ParsecLoader<T>::CompressMatrixHelper(configurations::Configurations &aConfigurations, std::unique_ptr<ExaGeoStatData<T>> &aData) {

    int max_rank = aConfigurations.GetMaxRank();
    int n = aConfigurations.GetProblemSize();
    int adaptive_decision = aConfigurations.GetAdaptiveDecision();
    int tol = aConfigurations.GetTolerance();
    int send_full_tile = 0;
    int auto_band = 0;
    int gpus = aConfigurations.GetGPUsNumbers();
    double upper_lower = EXAGEOSTAT_LOWER;
    int L = aConfigurations.GetDenseTileSize();
    int N = aConfigurations.GetProblemSize();
    int NT = (N % L == 0) ? (N/L) : (N/L + 1);

    SYNC_TIME_START();
    MatrixCompress((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), &ExaGeoStatHardware::GetHicmaParams()->norm_global, upper_lower, L, NT, max_rank, n,
                   adaptive_decision, tol, send_full_tile, auto_band, gpus, ExaGeoStatHardware::GetHicmaData(), ExaGeoStatHardware::GetParamsKernel());

    SYNC_TIME_PRINT(ExaGeoStatHardware::GetParsecMPIRank(), ("Matrix genneration Matrix norm: norm_global= %le\n", ExaGeoStatHardware::GetHicmaParams()->norm_global));
}

template<typename T>
void ParsecLoader<T>::ReadData(Configurations &aConfigurations, vector<T> &aMeasurementsMatrix, vector<T> &aXLocations,
                               vector<T> &aYLocations, vector<T> &aZLocations, const int &aP) {
}

template<typename T>
void ParsecLoader<T>::WriteData(const T &aMatrixPointer, const int &aProblemSize, const int &aP, std::string &aLoggerPath, Locations<T> &aLocations) {

}

template<typename T>
void ParsecLoader<T>::ReleaseInstance() {
    if (mpInstance != nullptr) {
        mpInstance = nullptr;
    }
}

template<typename T> ParsecLoader<T> *ParsecLoader<T>::mpInstance = nullptr;
