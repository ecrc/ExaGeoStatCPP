
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ParsecGenerator.cpp
 * @brief Implementation of the ParsecGenerator class
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2025-02-17
**/

#include <data-generators/concrete/runtime/ParsecGenerator.hpp>
#include <data-generators/LocationGenerator.hpp>
extern "C"{
    #include <runtime/parsec/jdf/JobDescriptionFormat.h>
}  

using namespace exageostat::common;
using namespace exageostat::results;
using namespace exageostat::configurations;
using namespace exageostat::generators::synthetic::parsec;

template<typename T>
ParsecGenerator<T> *ParsecGenerator<T>::GetInstance() {

    if (mpInstance == nullptr) {
        mpInstance = new ParsecGenerator<T>();
    }
    return mpInstance;
}

template<typename T>
std::unique_ptr<ExaGeoStatData<T>>
ParsecGenerator<T>::CreateSyntheticData(Configurations &aConfigurations, exageostat::kernels::Kernel<T> &aKernel) {

    // Initialize parameters
    int p = aKernel.GetVariablesNumber();
    int problem_size = aConfigurations.GetProblemSize() * p;
    auto data = std::make_unique<ExaGeoStatData<T>>(problem_size, aConfigurations.GetDimension());
    int HNB = aConfigurations.GetHNB();
    int dts = aConfigurations.GetDenseTileSize();
    int p_grid = aConfigurations.GetPGrid();
    int q_grid = aConfigurations.GetQGrid();
    int nodes = aConfigurations.GetCoresNumber();
    parsec_matrix_uplo_t uplo = PARSEC_MATRIX_LOWER;
    int rank = ExaGeoStatHardware::GetParsecMPIRank();
    int distance_matrix = aConfigurations.GetDistanceMetric();
    auto initial_theta = aConfigurations.GetInitialTheta();
    double sync_time_elapsed = 0.0;

    // Allocated new Locations object.
    auto *locations = new dataunits::Locations<T>(problem_size/p, aConfigurations.GetDimension());
    int parameters_number = aKernel.GetParametersNumbers();
    int seed = aConfigurations.GetSeed();
    int initial_seed[4] = {seed, seed, seed, 1};

    //normal random generation of e -- ei~problem_size(0, 1) to generate Z
    auto *randomN = new T[problem_size*100];
    LAPACKE_dlarnv(3, initial_seed, problem_size*100, (double *) randomN);

    // Set initial theta values.
    Configurations::InitTheta(aConfigurations.GetInitialTheta(), parameters_number);
    aConfigurations.SetInitialTheta(aConfigurations.GetInitialTheta());

    // Generate Locations phase
    LocationGenerator<T>::GenerateLocations(problem_size/p, aConfigurations.GetTimeSlot(), aConfigurations.GetDimension(),
                                            *locations);
    data->SetLocations(*locations);
    location l1;
    l1.x = (double *) locations->GetLocationX();
    l1.y = (double *) locations->GetLocationY();
    l1.z = (double *) locations->GetLocationZ();
    
    /* Data descritor declaration */
    data->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_Z);
    data->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_Z_COPY);
    data->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_DETERMINANT);
    data->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_PRODUCT);

    parsec_matrix_sym_block_cyclic_t *pDesc_C = new parsec_matrix_sym_block_cyclic_t;
    parsec_matrix_block_cyclic_t *pDesc_Z = data->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_Z).parsec_desc;
    parsec_matrix_block_cyclic_t *pDesc_Z_cpy = data->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_Z_COPY).parsec_desc;
    parsec_matrix_block_cyclic_t *pDesc_determinant = data->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_DETERMINANT).parsec_desc;
    parsec_matrix_block_cyclic_t *pDesc_product = data->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_PRODUCT).parsec_desc;

    auto *pParsec_context = (parsec_context_t *) ExaGeoStatHardware::GetParsecContext();

    parsec_matrix_sym_block_cyclic_init(pDesc_C, PARSEC_MATRIX_DOUBLE,
                                        rank, dts, dts, problem_size, problem_size, 0, 0,
                                        problem_size, problem_size, p_grid, nodes / p_grid, uplo);
    parsec_matrix_block_cyclic_init(pDesc_Z, PARSEC_MATRIX_DOUBLE, PARSEC_MATRIX_TILE,
                                    rank, dts, dts, problem_size, 1, 0, 0,
                                    problem_size, 1, 1, 1,
                                    1, 1, 0, 0);
    parsec_matrix_block_cyclic_init(pDesc_Z_cpy, PARSEC_MATRIX_DOUBLE, PARSEC_MATRIX_TILE,
                                    rank, dts, dts, problem_size, 1, 0, 0,
                                    problem_size, 1, 1, 1,
                                    1, 1, 0, 0);
    parsec_matrix_block_cyclic_init(pDesc_determinant, PARSEC_MATRIX_DOUBLE, PARSEC_MATRIX_TILE,
                                    rank, dts, dts, 1, 1, 0, 0,
                                    1, 1, 1, 1,
                                    1, 1, 0, 0);
    parsec_matrix_block_cyclic_init(pDesc_product, PARSEC_MATRIX_DOUBLE, PARSEC_MATRIX_TILE,
                                    rank, dts, dts, 1, 1, 0, 0,
                                    1, 1, 1, 1,
                                    1, 1, 0, 0);
    pDesc_C->mat = calloc((size_t) pDesc_C->super.nb_local_tiles *
                     (size_t) pDesc_C->super.bsiz,
                     (size_t) parsec_datadist_getsizeoftype(pDesc_C->super.mtype));
    parsec_data_collection_set_key((parsec_data_collection_t * )  pDesc_C, "pDesc_C");
    pDesc_Z->mat = calloc((size_t) pDesc_Z->super.nb_local_tiles *
                     (size_t) pDesc_Z->super.bsiz,
                     (size_t) parsec_datadist_getsizeoftype(pDesc_Z->super.mtype));
    parsec_data_collection_set_key((parsec_data_collection_t * )  pDesc_Z, "pDesc_Z");
    pDesc_Z_cpy->mat = calloc((size_t) pDesc_Z_cpy->super.nb_local_tiles *
                        (size_t) pDesc_Z_cpy->super.bsiz,
                        (size_t) parsec_datadist_getsizeoftype(pDesc_Z_cpy->super.mtype));
    parsec_data_collection_set_key((parsec_data_collection_t * )  pDesc_Z_cpy, "pDesc_Z_cpy");
    pDesc_determinant->mat = calloc((size_t) pDesc_determinant->super.nb_local_tiles *
                       (size_t) pDesc_determinant->super.bsiz,
                       (size_t) parsec_datadist_getsizeoftype(pDesc_determinant->super.mtype));
    parsec_data_collection_set_key((parsec_data_collection_t * )  pDesc_determinant, "pDesc_determinant");
    pDesc_product->mat = calloc((size_t) pDesc_product->super.nb_local_tiles *
                           (size_t) pDesc_product->super.bsiz,
                           (size_t) parsec_datadist_getsizeoftype(pDesc_product->super.mtype));
    parsec_data_collection_set_key((parsec_data_collection_t * )  pDesc_product, "pDesc_product");

    // Used to check results
    int info_potrf = 0;
    int info_trmm = 0;
    double sum = 0.0;

    /* Timer start */
    SYNC_TIME_START();

    ParsecDMatrixGeneration( pParsec_context, (parsec_tiled_matrix_t *)pDesc_C, &l1, &l1, &l1,
                               initial_theta.data(), distance_matrix, "univariate_matern_stationary", pDesc_C->super.lmt );
    SYNC_TIME_PRINT(rank, ("Matrix_generation_DOUBLE" "\tPxQ= %3d %-3d NB= %4d N= %7d HNB= %d : %lf\n",
                p_grid, nodes/p_grid, dts, problem_size, HNB, sync_time_elapsed));

    ParsecDMatrixSetDiagonal(pParsec_context, (parsec_tiled_matrix_t *)pDesc_C, aConfigurations.GetDiagonalAddition());
    /* Timer start */
    SYNC_TIME_START();
    ParsecDZGeneration( pParsec_context, (parsec_tiled_matrix_t *)pDesc_Z, (double *)&randomN[problem_size] ); //change i in C to 1
    SYNC_TIME_PRINT(rank, ("Z_vector_generation_DOUBLE" "\tPxQ= %3d %-3d NB= %4d N= %7d HNB= %d : %lf\n",
            p_grid, nodes/p_grid, dts, problem_size, HNB, sync_time_elapsed));

    sum = ParsecDMatrixSum(pParsec_context, (parsec_tiled_matrix_t *)(pDesc_C));
    fprintf( stderr, "\nin main C: sum after generation %.17g\n", sum );

    sum = ParsecDZSum(pParsec_context, (parsec_tiled_matrix_dc_t *)pDesc_Z);
    fprintf( stderr, "\nin main Z: sum after generation %.17g\n", sum );

    /* Free memory */
    parsec_data_free( pDesc_C->mat );
    parsec_tiled_matrix_dc_destroy( (parsec_tiled_matrix_t*)pDesc_C );

    auto hicma_data = ExaGeoStatHardware::GetHicmaData();
    {
        /* dcA0 */
        sym_two_dim_block_cyclic_init(&hicma_data->dcA0, (parsec_matrix_type_t ) matrix_RealDouble,
                rank, dts, dts, problem_size, problem_size, 0, 0,
                problem_size, problem_size, p_grid, nodes/p_grid, uplo);
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

    return data;
}

template<typename T>
void ParsecGenerator<T>::ReleaseInstance() {
    if (mpInstance != nullptr) {
        mpInstance = nullptr;
    }
}

template<typename T> ParsecGenerator<T> *ParsecGenerator<T>::mpInstance = nullptr;
