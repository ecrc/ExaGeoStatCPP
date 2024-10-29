
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file SyntheticGenerator.cpp
 * @brief Implementation of the SyntheticGenerator class
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-02-14
**/

#if !DEFAULT_RUNTIME
#include <runtime/parsec/ParsecHeader.h>

extern "C"{
#include <runtime/parsec/jdf/JobDescriptionFormat.h>
}

#endif

#include <data-loader/DataLoader.hpp>

using namespace std;

using namespace exageostat::dataLoader;
using namespace exageostat::common;
using namespace exageostat::results;

template<typename T>
std::unique_ptr<ExaGeoStatData<T>>
DataLoader<T>::CreateData(configurations::Configurations &aConfigurations, exageostat::kernels::Kernel<T> &aKernel) {

#if DEFAULT_RUNTIME
    // create vectors that will be populated with read data.
    vector<T> measurements_vector;
    vector<T> x_locations;
    vector<T> y_locations;
    vector<T> z_locations;

    aKernel.SetPValue(aConfigurations.GetTimeSlot());
    int p = aKernel.GetVariablesNumber();

    //Read the data out of the CSV file.
    this->ReadData(aConfigurations, measurements_vector, x_locations, y_locations, z_locations, p);

    //create data object
    auto data = std::make_unique<ExaGeoStatData<T>>(aConfigurations.GetProblemSize() / p,
                                                    aConfigurations.GetDimension());

    //Initialize the descriptors.
    auto linear_algebra_solver = linearAlgebra::LinearAlgebraFactory<T>::CreateLinearAlgebraSolver(EXACT_DENSE);

    linear_algebra_solver->InitiateDescriptors(aConfigurations, *data->GetDescriptorData(), p);
    linear_algebra_solver->ExaGeoStatLaSetTile(EXAGEOSTAT_UPPER_LOWER, 0, 0,
                                               data->GetDescriptorData()->GetDescriptor(CHAMELEON_DESCRIPTOR,
                                                                                        DESCRIPTOR_C).chameleon_desc);
    //populate data object with read data
    for (int i = 0; i < aConfigurations.GetProblemSize() / p; i++) {
        data->GetLocations()->GetLocationX()[i] = x_locations[i];
        data->GetLocations()->GetLocationY()[i] = y_locations[i];
        if (aConfigurations.GetDimension() != Dimension2D) {
            data->GetLocations()->GetLocationZ()[i] = z_locations[i];
        }
    }
    for (int i = 0; i < aConfigurations.GetProblemSize(); i++) {
        ((T *) data->GetDescriptorData()->GetDescriptor(CHAMELEON_DESCRIPTOR,
                                                        DESCRIPTOR_Z).chameleon_desc->mat)[i] = measurements_vector[i];
    }

    Results::GetInstance()->SetGeneratedLocationsNumber(aConfigurations.GetProblemSize() / p);
    Results::GetInstance()->SetIsLogger(aConfigurations.GetLogger());
    Results::GetInstance()->SetLoggerPath(aConfigurations.GetLoggerPath());
#else
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
    char filename[100];
    string files_directory_path = aConfigurations.GetDataPath();

    MB = L + 1;
    NB = L * 2;
	VERBOSE_PRINT(rank, verbose, ("Reading f_data\n"));
    data->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_F_DATA);
    parsec_matrix_block_cyclic_t *pF_data_desc = data->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_F_DATA).parsec_desc;
    ReadCSVToComplexTimeSlot((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), pF_data_desc, MB, NB, nodes, t, filename, rank, verbose, gpus);

    MB = 2*L-1;
    NB = L+1;
	VERBOSE_PRINT(rank, verbose, ("Reading Et1\n"));
    data->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_ET1);
    parsec_matrix_block_cyclic_t *pEt1_data_desc = data->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_ET1).parsec_desc;
    sprintf(filename, "%s/%s", files_directory_path.c_str(),"720_Et1.csv");
    ReadCSVComplex((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), pEt1_data_desc, MB, NB, nodes, t, filename, rank, verbose, gpus);

    MB = 2*L-1;
    NB = L-1;
	VERBOSE_PRINT(rank, verbose, ("Reading Et2\n"));
    data->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_ET2);
    parsec_matrix_block_cyclic_t *pEt2_data_desc = data->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_ET2).parsec_desc;
    sprintf(filename, "%s/%s", files_directory_path.c_str(),"720_Et2.csv");
    ReadCSVComplex((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), pEt2_data_desc, MB, NB, nodes, t, filename, rank, verbose, gpus);

    MB = 2*L;
    NB = 2*L-1;
	VERBOSE_PRINT(rank, verbose, ("Reading Ep\n"));
    data->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_EP);
    parsec_matrix_block_cyclic_t *pEp_data_desc = data->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_EP).parsec_desc;
    sprintf(filename, "%s/%s", files_directory_path.c_str(),"720_Ep.csv");
    ReadCSVComplex((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), pEp_data_desc, MB, NB, nodes, t, filename, rank, verbose, gpus);

    MB = (L*L+L)/2;
    NB = L;
	VERBOSE_PRINT(rank, verbose, ("Reading Slmn\n"));
    data->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_SLMN);
    parsec_matrix_block_cyclic_t *pSlum_data_desc = data->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_SLMN).parsec_desc;
    sprintf(filename, "%s/%s", files_directory_path.c_str(),"720_Slmn.csv");
    ReadCSVComplex((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), pSlum_data_desc, MB, NB, nodes, t, filename, rank, verbose, gpus);

    MB = L;
    NB = 2*L-1;
	VERBOSE_PRINT(rank, verbose, ("Reading Ie\n"));
    data->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_IE);
    parsec_matrix_block_cyclic_t *pIe_data_desc = data->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_IE).parsec_desc;
    sprintf(filename, "%s/%s", files_directory_path.c_str(),"720_Ie.csv");
    ReadCSVToComplex((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), pIe_data_desc, MB, NB, nodes, t, filename, rank, verbose, gpus);

    MB = L;
    NB = 2*L-1;
	VERBOSE_PRINT(rank, verbose, ("Reading Io\n"));
    data->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_IO);
    parsec_matrix_block_cyclic_t *pIo_data_desc = data->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_IO).parsec_desc;
    sprintf(filename, "%s/%s", files_directory_path.c_str(),"720_Io.csv");
    ReadCSVToComplex((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), pIo_data_desc, MB, NB, nodes, t, filename, rank, verbose, gpus);

    MB = L-1;
    NB = L+1;
	VERBOSE_PRINT(rank, verbose, ("Reading P\n"));
    data->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_P);
    parsec_matrix_block_cyclic_t *pP_data_desc = data->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_P).parsec_desc;
    sprintf(filename, "%s/%s", files_directory_path.c_str(),"720_P.csv");
    ReadCSVToComplex((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), pP_data_desc, MB, NB, nodes, t, filename, rank, verbose, gpus);

    MB = 2*L-1;
    NB = 2*L-1;
	VERBOSE_PRINT(rank, verbose, ("Reading D\n"));
    data->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_D);
    parsec_matrix_block_cyclic_t *pD_data_desc = data->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_D).parsec_desc;
    sprintf(filename, "%s/%s", files_directory_path.c_str(),"720_D.csv");
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
    sprintf(filename, "%s/%s", files_directory_path.c_str(),"720_flmERA.csv");
    ReadCSVTimeSlot((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), pFlmera_data_desc, MB, NB, nodes, t, filename, rank, verbose, gpus);

    // Backward
    if(aConfigurations.GetEnableInverse()){

        MB = L+1;
        NB = (L*L+L)/2;
        VERBOSE_PRINT(rank, verbose, ("Reading Zlm\n"));
        data->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_ZLM);
        parsec_matrix_block_cyclic_t *PZlm_data_desc = data->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_ZLM).parsec_desc;
        sprintf(filename, "%s/%s", files_directory_path.c_str(),"720_Zlm.csv");
        ReadCSV((parsec_context_t *) ExaGeoStatHardware::GetParsecContext(), PZlm_data_desc, MB, NB, nodes, t, filename, rank, verbose, gpus);

        MB = 2*L-1;
        NB = 2*L;
        VERBOSE_PRINT(rank, verbose, ("Reading SC\n"));
        data->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_SC);
        parsec_matrix_block_cyclic_t *pSc_data_desc = data->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_SC).parsec_desc;
        sprintf(filename, "%s/%s", files_directory_path.c_str(),"720_SC.csv");
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
	VERBOSE_PRINT(rank, verbose, ("Alocating FLMT\n"));
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
        sprintf(filename, "%s/%s", files_directory_path.c_str(),"ts_test.csv");
        ReadCSVFileHelper(filename, pFileContent, ts_test_M, ts_test_N);
    }

#endif
    return data;
}


template<typename T>
int DataLoader<T>::ReadCSVFileHelper(const char* apFilename, double *apFileContent, int aM, int aN) {

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