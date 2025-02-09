
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file SyntheticGenerator.cpp
 * @brief Implementation of the SyntheticGenerator class
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2024-02-04
**/

#include <data-generators/concrete/SyntheticGenerator.hpp>
#include <data-generators/LocationGenerator.hpp>
#if !DEFAULT_RUNTIME
#include <data-loader/concrete/ParsecLoader.hpp>
extern "C"{
#include <runtime/parsec/jdf/JobDescriptionFormat.h>
}
#else
#include <data-loader/concrete/CSVLoader.hpp>
#endif
//TODO: we need to make WriteData a function outside the csv, So it can be used whatever the runtime is.
// currently, it has an implementation for the CSVLoader and an empty body for the parsec loader
using namespace exageostat::generators::synthetic;
using namespace exageostat::common;
using namespace exageostat::configurations;
using namespace exageostat::results;

template<typename T>
SyntheticGenerator<T> *SyntheticGenerator<T>::GetInstance() {

    if (mpInstance == nullptr) {
        mpInstance = new SyntheticGenerator<T>();
    }
    return mpInstance;
}

template<typename T>
std::unique_ptr<ExaGeoStatData<T>>
SyntheticGenerator<T>::CreateData(Configurations &aConfigurations,
                                  exageostat::kernels::Kernel<T> &aKernel) {
#if DEFAULT_RUNTIME
    int n = aConfigurations.GetProblemSize() * aConfigurations.GetTimeSlot();
    auto data = std::make_unique<ExaGeoStatData<T>>(n, aConfigurations.GetDimension());

    // Allocated new Locations object.
    auto *locations = new dataunits::Locations<T>(n, aConfigurations.GetDimension());
    int parameters_number = aKernel.GetParametersNumbers();

    // Set initial theta values.
    Configurations::InitTheta(aConfigurations.GetInitialTheta(), parameters_number);
    aConfigurations.SetInitialTheta(aConfigurations.GetInitialTheta());

    // Generate Locations phase
    LocationGenerator<T>::GenerateLocations(n, aConfigurations.GetTimeSlot(), aConfigurations.GetDimension(),
                                            *locations);
    data->SetLocations(*locations);

    // TODO: May need to get refactored to avoid the if/else guards

    // Generate Descriptors phase
    auto linear_algebra_solver = linearAlgebra::LinearAlgebraFactory<T>::CreateLinearAlgebraSolver(EXACT_DENSE);
    linear_algebra_solver->GenerateSyntheticData(aConfigurations, data, aKernel);

    if (aConfigurations.GetLogger()) {
        VERBOSE("Writing generated data to the disk (Synthetic Dataset Generation Phase) .....")
#ifdef USE_MPI
        auto pMatrix = new T[aConfigurations.GetProblemSize()];
        std::string path = aConfigurations.GetLoggerPath();

        CHAMELEON_Desc2Lap(ChamUpperLower, data->GetDescriptorData()->GetDescriptor(CHAMELEON_DESCRIPTOR,
                                                                                    DESCRIPTOR_Z).chameleon_desc,
                           pMatrix, aConfigurations.GetProblemSize());
        if (helpers::CommunicatorMPI::GetInstance()->GetRank() == 0) {
            dataLoader::csv::CSVLoader<T>::GetInstance()->WriteData(
                    *((T *) data->GetDescriptorData()->GetDescriptor(CHAMELEON_DESCRIPTOR,
                                                                     DESCRIPTOR_Z).chameleon_desc->mat),
                    aConfigurations.GetProblemSize(), parameters_number, path,
                    *data->GetLocations());
        }
        delete[] pMatrix;
#else
        std::string path = aConfigurations.GetLoggerPath();
        dataLoader::csv::CSVLoader<T>::GetInstance()->WriteData(*((T *) data->GetDescriptorData()->GetDescriptor(CHAMELEON_DESCRIPTOR,
                                                                                           DESCRIPTOR_Z).chameleon_desc->mat),
                                          aConfigurations.GetProblemSize(), aKernel.GetVariablesNumber(), path,
                                          *data->GetLocations());
#endif
    }
#else


    int p = aKernel.GetVariablesNumber();
    int N = aConfigurations.GetProblemSize() * p;
    auto data = std::make_unique<ExaGeoStatData<T>>(N, aConfigurations.GetDimension());
    int HNB = aConfigurations.GetHNB();
    int dts = aConfigurations.GetDenseTileSize();
    int p_grid = aConfigurations.GetPGrid();
    int q_grid = aConfigurations.GetQGrid();
    int nodes = aConfigurations.GetCoresNumber();
    parsec_matrix_uplo_t uplo = PARSEC_MATRIX_LOWER;
    int rank = ExaGeoStatHardware::GetParsecMPIRank();
    int distance_matrix = aConfigurations.GetDistanceMetric();
    auto initial_theta = aConfigurations.GetInitialTheta();

    // Allocated new Locations object.
    auto *locations = new dataunits::Locations<T>(N/p, aConfigurations.GetDimension());
    int parameters_number = aKernel.GetParametersNumbers();
    int seed = aConfigurations.GetSeed();
    int initial_seed[4] = {seed, seed, seed, 1};


    //normal random generation of e -- ei~N(0, 1) to generate Z
    auto *randomN = new T[N*100];
    LAPACKE_dlarnv(3, initial_seed, N*100, (double *) randomN);

    // Set initial theta values.
    Configurations::InitTheta(aConfigurations.GetInitialTheta(), parameters_number);
    aConfigurations.SetInitialTheta(aConfigurations.GetInitialTheta());

    // Generate Locations phase
    LocationGenerator<T>::GenerateLocations(N/p, aConfigurations.GetTimeSlot(), aConfigurations.GetDimension(),
                                            *locations);
    data->SetLocations(*locations);
    location l1;
    l1.x = (double *) locations->GetLocationX();
    l1.y = (double *) locations->GetLocationY();
    l1.z = (double *) locations->GetLocationZ();
    /* Data descritor declaration */
    parsec_matrix_sym_block_cyclic_t dcC;
    parsec_matrix_block_cyclic_t dcZ;
    parsec_matrix_block_cyclic_t dcZcpy;
    parsec_matrix_block_cyclic_t dcdet;
    parsec_matrix_block_cyclic_t dcproduct;

    auto *parsec = (parsec_context_t *) ExaGeoStatHardware::GetParsecContext();

    double sync_time_elapsed = 0.0;

    parsec_matrix_sym_block_cyclic_init(&dcC, PARSEC_MATRIX_DOUBLE,
                                        rank, dts, dts, N, N, 0, 0,
                                        N, N, p_grid, nodes / p_grid, uplo);
    parsec_matrix_block_cyclic_init(&dcZ, PARSEC_MATRIX_DOUBLE, PARSEC_MATRIX_TILE,
                                    rank, dts, dts, N, 1, 0, 0,
                                    N, 1, 1, 1,
                                    1, 1, 0, 0);
    parsec_matrix_block_cyclic_init(&dcZcpy, PARSEC_MATRIX_DOUBLE, PARSEC_MATRIX_TILE,
                                    rank, dts, dts, N, 1, 0, 0,
                                    N, 1, 1, 1,
                                    1, 1, 0, 0);
    parsec_matrix_block_cyclic_init(&dcdet, PARSEC_MATRIX_DOUBLE, PARSEC_MATRIX_TILE,
                                    rank, dts, dts, 1, 1, 0, 0,
                                    1, 1, 1, 1,
                                    1, 1, 0, 0);
    parsec_matrix_block_cyclic_init(&dcproduct, PARSEC_MATRIX_DOUBLE, PARSEC_MATRIX_TILE,
                                    rank, dts, dts, 1, 1, 0, 0,
                                    1, 1, 1, 1,
                                    1, 1, 0, 0);
    dcC.mat = calloc((size_t) dcC.super.nb_local_tiles *
                     (size_t) dcC.super.bsiz,
                     (size_t) parsec_datadist_getsizeoftype(dcC.super.mtype));
    parsec_data_collection_set_key((parsec_data_collection_t * ) & dcC, "dcC");
    dcZ.mat = calloc((size_t) dcZ.super.nb_local_tiles *
                     (size_t) dcZ.super.bsiz,
                     (size_t) parsec_datadist_getsizeoftype(dcZ.super.mtype));
    parsec_data_collection_set_key((parsec_data_collection_t * ) & dcZ, "dcZ");
    dcZcpy.mat = calloc((size_t) dcZcpy.super.nb_local_tiles *
                        (size_t) dcZcpy.super.bsiz,
                        (size_t) parsec_datadist_getsizeoftype(dcZcpy.super.mtype));
    parsec_data_collection_set_key((parsec_data_collection_t * ) & dcZcpy, "dcZcpy");
    dcdet.mat = calloc((size_t) dcdet.super.nb_local_tiles *
                       (size_t) dcdet.super.bsiz,
                       (size_t) parsec_datadist_getsizeoftype(dcdet.super.mtype));
    parsec_data_collection_set_key((parsec_data_collection_t * ) & dcdet, "dcdet");
    dcproduct.mat = calloc((size_t) dcproduct.super.nb_local_tiles *
                           (size_t) dcproduct.super.bsiz,
                           (size_t) parsec_datadist_getsizeoftype(dcproduct.super.mtype));
    parsec_data_collection_set_key((parsec_data_collection_t * ) & dcproduct, "dcproduct");

    printf("\n");
    for (int i = 0; i < N; i++) {
        printf("x: %f, y: %f\n", l1.x[i], l1.y[i]);
    }

    // Used to check results
    int info_potrf = 0;
    int info_trmm = 0;
    double sum = 0.0;

    /* Timer start */
    SYNC_TIME_START();

    ParsecDMatrixGeneration( parsec, (parsec_tiled_matrix_t *)&dcC, &l1, &l1, &l1,
                               initial_theta.data(), distance_matrix, "univariate_matern_stationary", dcC.super.lmt );
    SYNC_TIME_PRINT(rank, ("Matrix_generation_DOUBLE" "\tPxQ= %3d %-3d NB= %4d N= %7d HNB= %d : %lf\n",
                p_grid, nodes/p_grid, dts, N, HNB, sync_time_elapsed));

    ParsecDMatrixSetDiagonal(parsec, (parsec_tiled_matrix_t *)&dcC, aConfigurations.GetDiagonalAddition());
    /* Timer start */
    SYNC_TIME_START();
    // TODO: 99 for sure is wrong!!
    int i = 99;
    printf("i? %d, %d\n", i, i*N);
    ParsecDZGeneration( parsec, (parsec_tiled_matrix_t *)&dcZ, (double *)&randomN[i*N] );
    SYNC_TIME_PRINT(rank, ("Z_vector_generation_DOUBLE" "\tPxQ= %3d %-3d NB= %4d N= %7d HNB= %d : %lf\n",
            p_grid, nodes/p_grid, dts, N, HNB, sync_time_elapsed));

    sum = ParsecDMatrixSum(parsec, (parsec_tiled_matrix_t *)&(dcC));
    fprintf( stderr, "\nin main C: sum after generation %.17g\n", sum );

    sum = ParsecDZSum(parsec, (parsec_tiled_matrix_dc_t *)&dcZ);
    fprintf( stderr, "\nin main Z: sum after generation %.17g\n", sum );

            /* Free memory */
    parsec_data_free( dcC.mat );
    parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_t*)&dcC );

    auto hicma_data = ExaGeoStatHardware::GetHicmaData();
    {
        /* dcA0 */
        sym_two_dim_block_cyclic_init(&hicma_data->dcA0, (parsec_matrix_type_t ) matrix_RealDouble,
                rank, dts, dts, N, N, 0, 0,
                N, N, p_grid, nodes/p_grid, uplo);
        hicma_data->dcA0.mat = parsec_data_allocate((size_t)hicma_data->dcA0.super.nb_local_tiles *
                (size_t)hicma_data->dcA0.super.bsiz *
                (size_t)parsec_datadist_getsizeoftype(hicma_data->dcA0.super.mtype));
        parsec_data_collection_set_key((parsec_data_collection_t*)&hicma_data.dcA0, "dcA0");
    }
    int ret = 0;

            /* Check results */
        if( 0 == rank && info_potrf != 0 ) {
            fprintf(stderr, "-- Cholesky factorization is suspicious (info_potrf = %d) ! \n", info_potrf);
            ret |= 1;
            //exit(1);
        }

        if( 0 == rank && info_trmm != 0 ) {
            fprintf(stderr, "-- TRMM is suspicious (info_potrf = %d) ! \n", info_trmm);
            ret |= 1;
        }

        if( 0 == rank ) {
            fprintf(stderr, "Done Z Vector Generation Phase.\n");
            fprintf(stderr, "************************************************************\n");
        }

    Results::GetInstance()->SetGeneratedLocationsNumber(aConfigurations.GetProblemSize());
    Results::GetInstance()->SetIsLogger(aConfigurations.GetLogger());
    Results::GetInstance()->SetLoggerPath(aConfigurations.GetLoggerPath());
#endif

    VERBOSE("Done.")
    return data;
}

template<typename T>
void SyntheticGenerator<T>::ReleaseInstance() {
    if (mpInstance != nullptr) {
        mpInstance = nullptr;
    }
}

template<typename T> SyntheticGenerator<T> *SyntheticGenerator<T>::mpInstance = nullptr;
