
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
 * @author Rabab Alomairy
 * @date 2025-11-13
**/

#include <nlopt.hpp>

#include <runtime-solver/concrete/ParsecRuntimeSolver.hpp>
#include <data-analyzer/DataAnalyzer.hpp>
#include <data-transformer/DataTransformer.hpp>
#include <results/Results.hpp>

#include <runtime/parsec/ParsecHeader.h>
extern "C"{
#include <runtime/parsec/jdf/JobDescriptionFormat.h>
}

#include <data-units/ModelingDataHolders.hpp>
#include <utilities/Logger.hpp>

using namespace std;
using namespace nlopt;

using namespace exageostat::common;
using namespace exageostat::configurations;
using namespace exageostat::analyzer;
using namespace exageostat::dataunits;
using namespace exageostat::runtimesolver;
using namespace exageostat::results;

template<typename T>
void ParsecRuntimeSolver<T>::ExaGeoStatSYRK(unique_ptr<ExaGeoStatData<T>> &aData){
    auto* pContext = (parsec_context_t *) ExaGeoStatHardware::GetParsecContext();
    auto* pDesc_A = (parsec_tiled_matrix_t *) aData->GetDescriptorData()->GetDescriptor(DescriptorType::PARSEC_DESCRIPTOR, DescriptorName::DESCRIPTOR_A).parsec_desc;

    SYNC_TIME_START();
    dplasma_dsyrk(pContext, dplasmaLower, dplasmaNoTrans, 1.0, pDesc_A, 0.0, (parsec_tiled_matrix_t *) &ExaGeoStatHardware::GetHicmaData()->dcA);
    SYNC_TIME_PRINT(ExaGeoStatHardware::GetParsecMPIRank(),("SYRK\n"));
}

template<typename T>
void ParsecRuntimeSolver<T>::ExaGeoStatTLRCholesky(unique_ptr<ExaGeoStatData<T>> &aData){

    auto *pParams = ExaGeoStatHardware::GetHicmaParams();
    auto *pHicma_data = ExaGeoStatHardware::GetHicmaData();
    auto *pAnalysis = ExaGeoStatHardware::GetAnalysis();
    auto *pContext = (parsec_context_t *) ExaGeoStatHardware::GetParsecContext();
    for( int i= 0; i < pParams->nruns; i++ ) {
        hicma_parsec_potrf(pContext, pHicma_data, pParams, pAnalysis);
    }
}

template<typename T>
double ParsecRuntimeSolver<T>::ExaGeoStatNorm(Configurations &aConfigurations, unique_ptr<ExaGeoStatData<T>> &aData){

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
double ParsecRuntimeSolver<T>::CalculateMSE(Configurations &aConfigurations, unique_ptr<ExaGeoStatData<T>> &aData) {

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
T ParsecRuntimeSolver<T>::ModelingOperations(unique_ptr<ExaGeoStatData<T>> &aData, Configurations &aConfigurations,
                                              T *apMeasurementsMatrix, const kernels::Kernel<T> &aKernel) {

    // Route based on IsClimateEmulator flag
    if (aConfigurations.GetIsClimateEmulator()) {
        // ============ CLIMATE EMULATOR PATH ============
        // Direct matrix operations (no MLE optimization)
        // Uses SYRK + Cholesky decomposition for climate data
        
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
        
        
    } else {
        // ============ GENERAL MLE OPTIMIZATION PATH ============
        // Standard geospatial modeling with MLE parameter estimation
        
        int parameters_number = aKernel.GetParametersNumbers();
        int max_number_of_iterations = aConfigurations.GetMaxMleIterations();
        // Setting struct of data to pass to the modeling.
        auto modeling_data = new mModelingData(aData, aConfigurations, *apMeasurementsMatrix, aKernel);
        // Create nlopt
        double opt_f;
        opt optimizing_function(nlopt::LN_BOBYQA, parameters_number);
        // Initialize problem's bound.
        optimizing_function.set_lower_bounds(aConfigurations.GetLowerBounds());
        optimizing_function.set_upper_bounds(aConfigurations.GetUpperBounds());
        optimizing_function.set_ftol_abs(aConfigurations.GetTolerance());
        // Set max iterations value.
        optimizing_function.set_maxeval(max_number_of_iterations);
        optimizing_function.set_max_objective(DataModelingAPI, (void *) modeling_data);
        
        // Optimize mle using nlopt.
        optimizing_function.optimize(aConfigurations.GetStartingTheta(), opt_f);
        aConfigurations.SetEstimatedTheta(aConfigurations.GetStartingTheta());
        
        auto theta = aConfigurations.GetStartingTheta();
        
        LOGGER("--> Final Theta Values (", true)
        for (int i = 0; i < parameters_number; i++) {
            LOGGER_PRECISION(theta[i])
            if (i != parameters_number - 1) {
                LOGGER_PRECISION(", ")
            }
        }
        LOGGER_PRECISION(")")
        LOGGER("")
        
        // delete modeling_data;
        return optimizing_function.last_optimum_value();
    }

}

template<typename T>
double ParsecRuntimeSolver<T>::DataModelingAPI(const vector<double> &aTheta, vector<double> &aGrad, void *apInfo) {

    auto config = ((mModelingData<T> *) apInfo)->mpConfiguration;
    auto data = ((mModelingData<T> *) apInfo)->mpData;
    auto measurements = ((mModelingData<T> *) apInfo)->mpMeasurementsMatrix;
    auto kernel = ((mModelingData<T> *) apInfo)->mpKernel;
   
    return ParsecDmleTile(*config,aTheta.data(),aGrad,*data,*kernel) ; 
}

template<typename T>
double ParsecRuntimeSolver<T>::ParsecDmleTile(Configurations &aConfigurations,const double *theta, vector<double> & grad,
    unique_ptr<ExaGeoStatData<T>> &aData,const kernels::Kernel<T> &aKernel) {

    // Initialization of timing and performance variables.
    T loglik = 0.0, logdet = 0.0, time_facto = 0.0, time_solve = 0.0, time_mle = 0.0, variance = 1;
    T logdet_calculate = 0.0, matrix_gen_time = 0.0, dzcpy_time = 0.0;
    int N, NRHS, success;
    T flops = 0.0;
    double accumulated_executed_time, accumulated_flops;
    
    int p = aKernel.GetVariablesNumber();
    int problem_size = aConfigurations.GetProblemSize() * p;
    int HNB = aConfigurations.GetHNB();
    int dts = aConfigurations.GetDenseTileSize();
    int p_grid = aConfigurations.GetPGrid();
    int q_grid = aConfigurations.GetQGrid();
    int nodes = aConfigurations.GetCoresNumber();
    parsec_matrix_uplo_t uplo = PARSEC_MATRIX_LOWER;
    int rank = ExaGeoStatHardware::GetParsecMPIRank();
    int distance_matrix = aConfigurations.GetDistanceMetric();
    auto initial_theta = aConfigurations.GetInitialTheta();

    // Get the MLE data structure.
    auto kernel_name = aConfigurations.GetKernelName();
    int num_params = aKernel.GetParametersNumbers();
    
    T *determinant = aData->GetDescriptorData()->GetDescriptorMatrix(PARSEC_DESCRIPTOR, DESCRIPTOR_DETERMINANT);
    *determinant = 0;
    T *product = aData->GetDescriptorData()->GetDescriptorMatrix(PARSEC_DESCRIPTOR, DESCRIPTOR_PRODUCT);
    *product = 0;
    double sum = 0.0;
   
    hicma_parsec_params_t* hicma_params = ExaGeoStatHardware::GetHicmaParams();
    hicma_parsec_data_t* hicma_data = ExaGeoStatHardware::GetHicmaData();
    starsh_params_t* params_kernel = ExaGeoStatHardware::GetParamsKernel();
    hicma_parsec_matrix_analysis_t* analysis = ExaGeoStatHardware::GetAnalysis();
    hicma_params->norm_global = 0.0;
    double threshold = pow(10.0, -1.0 * 16);

#if SWITCH_TO_DP
    double threshold_to_dp = pow(10.0, -1.0 * data->opt_tol_mix);
#endif

auto *parsec = (parsec_context_t *) ExaGeoStatHardware::GetParsecContext();

// Temporarily generate data: l1 struct is not part of our Locations class
    auto locations = aData->GetLocations();
    location l1;
    l1.x = (double *) locations->GetLocationX();
    l1.y = (double *) locations->GetLocationY();
    l1.z = (double *) locations->GetLocationZ();

// end of temporarily data generation

    auto* descC =(parsec_tiled_matrix_t *) &hicma_data->dcA;
    auto* descZ = &((aData->GetDescriptorData()->GetDescriptor(DescriptorType::PARSEC_DESCRIPTOR, DescriptorName::DESCRIPTOR_Z).parsec_desc)->super);
    auto* descZcpy = &((aData->GetDescriptorData()->GetDescriptor(DescriptorType::PARSEC_DESCRIPTOR, DescriptorName::DESCRIPTOR_Z_COPY).parsec_desc)->super);
    auto* descdet = &((aData->GetDescriptorData()->GetDescriptor(DescriptorType::PARSEC_DESCRIPTOR, DescriptorName::DESCRIPTOR_DETERMINANT).parsec_desc)->super);
    auto* descproduct = &((aData->GetDescriptorData()->GetDescriptor(DescriptorType::PARSEC_DESCRIPTOR, DescriptorName::DESCRIPTOR_PRODUCT).parsec_desc)->super);

    N = descC->m;
    NRHS = descZ->n;

    string recovery_file = aConfigurations.GetRecoveryFile();
    int iter_count = aData->GetMleIterations();

    START_TIMING(dzcpy_time);
    if (iter_count == 0) {
        VERBOSE("Save copy of the original Z vector...");
        dplasma_dlacpy(parsec, PlasmaUpperLower, descZ, descZcpy);
        VERBOSE(" Done.\n");
    }
    STOP_TIMING(dzcpy_time);

   { START_TIMING(dzcpy_time);
    if(iter_count > 0) {
        VERBOSE("Re-store the original Z vector...");
        dplasma_dlacpy(parsec, PlasmaUpperLower, descZcpy, descZ);
        VERBOSE(" Done.\n");
    }
    STOP_TIMING(dzcpy_time);
   }

#if DEBUG_INFO
        VERBOSE("Get sum of C...");
        sum = ParsecDMatrixSum(parsec, descC);
        if (rank == 0)
            fprintf(stderr, "\nC: sum before generation %.17g\n", sum);
        VERBOSE("Get sum of Z...");
        sum = ParsecDZSum(parsec, descZ);
        if (rank == 0)
            fprintf(stderr, "\nZ: sum before generation %.17g\n", sum);
#endif

#if defined(EXAGEOSTAT_USE_CUDA)
        parsec_devices_release_memory();
        parsec_devices_reset_load(parsec);
#endif

        // Initialize decisions for each tile.
        hicma_parsec_decision_init(hicma_params);

        VERBOSE("Generate New Covariance Matrix...");
        START_TIMING(matrix_gen_time);
        SYNC_TIME_START();
        STARSH_ssdata* sdata = (STARSH_ssdata*)(params_kernel->data);
        sdata->sigma = theta[0];
        sdata->beta  = theta[1];
        sdata->nu    = theta[2];
        sdata->beta_time = theta[3];
        sdata->nu_time   = theta[4];
        sdata->nonsep_param = theta[5];
        printf("**** thetas: %f %f %f %f %f %f\n",theta[0],theta[1],theta[2],theta[3],theta[4],theta[5]);
        // Generate the covariance matrix.
        hicma_parsec_matrix_generation(parsec, hicma_data, hicma_params, params_kernel);
        SYNC_TIME_PRINT(hicma_params->rank, ("Matrix compression\tband_size_norm= %3d norm_global= %le\n",
                                             hicma_params->band_size_norm, hicma_params->norm_global));
        hicma_params->time_starsh = sync_time_elapsed;
        STOP_TIMING(matrix_gen_time);
        VERBOSE(" Done.\n");

#if DEBUG_INFO
        SYNC_TIME_START();
        hicma_parsec_rank_check(parsec, (parsec_tiled_matrix_t*)(&hicma_data->dcAr),
                                hicma_params->band_size_dense);
        SYNC_TIME_PRINT(hicma_params->rank, ("Check Rank after compression\n"));
#endif

        // Gathering rank info.
        SYNC_TIME_START();
        hicma_parsec_rank_stat(parsec, "init_rank_tile", (parsec_tiled_matrix_t*)(&hicma_data->dcAr),
                                 hicma_params, &hicma_params->iminrk, &hicma_params->imaxrk, &hicma_params->iavgrk);
        SYNC_TIME_PRINT(hicma_params->rank, ("Gather rank after matrix compression\tmin= %d max= %d avg= %lf\n",
                                             hicma_params->iminrk, hicma_params->imaxrk, hicma_params->iavgrk));

        // Band size auto-tuning.
        SYNC_TIME_START();
        hicma_parsec_band_size_dense_auto_tuning(parsec, hicma_data, hicma_params, params_kernel);
        SYNC_TIME_PRINT(hicma_params->rank, ("Band_size auto-tuning\tband_size_dense= %d band_size_dist= %d band_size_auto_tuning_termination= %lf\n",
                                             hicma_params->band_size_dense, hicma_params->band_size_dist,
                                             hicma_params->band_size_auto_tuning_termination));
        hicma_params->time_opt_band = sync_time_elapsed;
        if (hicma_params->auto_band && hicma_params->band_size_dense > 1) {
            SYNC_TIME_START();
            hicma_parsec_rank_stat(parsec, "init_rank_tile_auto_band",
                                     reinterpret_cast<parsec_tiled_matrix_t*>(&hicma_data->dcAr),
                                     hicma_params, &hicma_params->iminrk_auto_band,
                                     &hicma_params->imaxrk_auto_band, &hicma_params->iavgrk_auto_band);
            SYNC_TIME_PRINT(hicma_params->rank, ("Gather rank after band_size auto-tuning\tmin= %d max= %d avg= %lf\n",
                                                 hicma_params->iminrk_auto_band, hicma_params->imaxrk_auto_band,
                                                 hicma_params->iavgrk_auto_band));
        }

        // Memory calculation.
        SYNC_TIME_START();
        hicma_parsec_memory_calculation(parsec, hicma_params);
        SYNC_TIME_PRINT(hicma_params->rank, ("Memory for matrix allocation\tactual= %lf GB maxrank= %lf GB\n",
                                             hicma_params->memory_per_node, hicma_params->memory_per_node_maxrank));

        // Sparse analysis.
        SYNC_TIME_START();
        hicma_parsec_sparse_analysis(parsec, hicma_data, hicma_params, analysis);
        SYNC_TIME_PRINT(hicma_params->rank, ("Matrix analysis\tsparse= %d band_size_dist= %d WORKLOAD_BALANCE= %d\n",
                                             hicma_params->sparse, hicma_params->band_size_dist, WORKLOAD_BALANCE));
        hicma_params->time_analysis = sync_time_elapsed;

        if (hicma_params->check) {
            SYNC_TIME_START();
            hicma_parsec_check_compression(parsec, hicma_data, hicma_params, params_kernel, analysis);
            SYNC_TIME_PRINT(hicma_params->rank, ("Check compress and run dense dpotrf\n"));
        }

#if DEBUG_INFO
        SYNC_TIME_START();
        hicma_parsec_rank_check(parsec, reinterpret_cast<parsec_tiled_matrix_t*>(&hicma_data->dcAr),
                                hicma_params->band_size_dense);
        SYNC_TIME_PRINT(hicma_params->rank, ("Before HiCMA Check Rank\n"));
#endif

#if PRINT_RANK
        hicma_parsec_process_id_print(hicma_data, hicma_params);
#endif

#if COUNT_VALUE_HALF 
        VERBOSE("Count values exceeding half-precision...");
        if (data->iter_count == 0)
            parsec_dmatrix_ratio_error_half(parsec, descC, theta, iter_count);
        double count_value_time = 0.0;
        START_TIMING(count_value_time);
        parsec_dmatrix_ratio_exceed_half(parsec, descC, theta, iter_count);
        STOP_TIMING(count_value_time);
        VERBOSE(" Done.\n");
#if defined(EXAGEOSTAT_USE_CUDA)
        parsec_devices_release_memory();
        parsec_devices_reset_load(parsec);
#endif
#endif

        VERBOSE("Cholesky factorization of Sigma...");
        START_TIMING(time_facto);
#if DEBUG_INFO
        sum = ParsecDMatrixSum(parsec, (parsec_tiled_matrix_t *)(dcC));
        if (rank == 0)
            fprintf(stderr, "\nC: sum before cholesky %.17g\n", sum);
#endif

    /* HiCMA Cholesky */
        SYNC_TIME_START();
        hicma_params->info = hicma_parsec_potrf(parsec, hicma_data, hicma_params, analysis);
        SYNC_TIME_PRINT(hicma_params->rank, ("hicma_parsec_cholesky here\tband_size_dense= %d band_size_dist= %d lookahead= %d kind_of_problem= %d HNB= %d PxQ= %3d %-3d nb_gpus= %d NB= %4d N= %7d kind_of_cholesky= %d sparse= %d nb_dense_dp= %.0lf nb_dense_sp= %.0lf nb_dense_hp= %.0lf nb_low_rank_dp= %.0lf nb_low_rank_sp= %.0lf: %14f gflops\n",
        hicma_params->band_size_dense, hicma_params->band_size_dist,
        hicma_params->lookahead, hicma_params->kind_of_problem, hicma_params->HNB,
        hicma_params->P, hicma_params->Q, hicma_params->gpus, hicma_params->NB,
        hicma_params->N, hicma_params->kind_of_cholesky, hicma_params->sparse,
        hicma_params->nb_dense_dp, hicma_params->nb_dense_sp,
        hicma_params->nb_dense_hp, hicma_params->nb_low_rank_dp, hicma_params->nb_low_rank_sp,
        hicma_params->gflops = (hicma_params->flops / 1e9) / sync_time_elapsed));
        hicma_params->time_hicma = sync_time_elapsed;
#if DEBUG_INFO
        dplasma_dprint(parsec, hicma_params->uplo, (parsec_tiled_matrix_t*)(&hicma_data->dcA));
#endif

    if (rank == 0 && hicma_params->info != 0) {
        fprintf(stderr, "-- Factorization is suspicious (info = %d) ! \n", hicma_params->info);
    }
    
    auto NT = descC->lmt;
    auto* rank_array = hicma_params->rank_array;
     analysis->trsm_initial = (uint16_t **)malloc( NT * sizeof(uint16_t *) );
    for(int i = 0; i < NT-1; i++)
        analysis->trsm_initial[i] = (uint16_t *)malloc( (NT-i+1) * sizeof(uint16_t) );
    analysis->trsm_num_initial = (uint16_t *)calloc( NT, sizeof(uint16_t) );

    /* Allocate memory */
    /* Final rank ditribution */
    analysis->initial_rank = (uint16_t *)calloc( NT * NT, sizeof(uint16_t) );
    analysis->final_rank = (uint8_t *)calloc( NT * NT, sizeof(uint8_t) );
    for(int j = 0; j < NT-1; j++) {
        int num_tmp = 0;
        for(int i = j+1; i < NT; i++) {
            if( rank_array[j*NT+i] > 0 ) {
                analysis->initial_rank[j*NT+i] = (uint16_t)rank_array[j*NT+i];
                analysis->final_rank[j*NT+i] = (uint8_t)1;
                analysis->initial_density += 1.0;
                analysis->trsm_initial[j][num_tmp] = (uint16_t)i;
                num_tmp += 1;
            }
        }
        analysis->trsm_num_initial[j] = (uint16_t)num_tmp;
    }
    analysis->initial_density = analysis->initial_density / NT / NT * 2;
    hicma_parsec_matrix_uncompress(parsec, uplo,(parsec_tiled_matrix_t*)(&hicma_data->dcA0),
    (parsec_tiled_matrix_t*)(&hicma_data->dcA), (parsec_tiled_matrix_t*)(&hicma_data->dcAr),
    analysis, hicma_params->band_size_dense, hicma_params->maxrank, &hicma_params->info);
    if (hicma_params->info != 0) {
        if (rank == 0)
            fprintf(stderr, "\nCheck Uncompression : %d\n", hicma_params->info);
        exit(1);
    }

#if DEBUG_INFO
for(int i = 0; i< hicma_data->dcA0.grid.rows; i++){

    for(int j = 0; j< hicma_data->dcA0.grid.cols; j++){
        printf("%f ", ((double *) hicma_data->dcA0.mat)[i * hicma_data->dcA0.grid.cols +j]);
    }
    printf("\n");
}
        sum = ParsecDMatrixSum(parsec, (parsec_tiled_matrix_t *)(&hicma_data->dcA0));
        if (rank == 0)
            fprintf(stderr, "\nC: sum after cholesky %.17g\n", sum);
#endif

        VERBOSE("Calculating the log determinant ...");
        START_TIMING(logdet_calculate);
        ParsecDMatrixDet(parsec, (parsec_tiled_matrix_t*)(&hicma_data->dcA0), descdet);
        int root = descdet->super.rank_of(&descdet->super, 0, 0);
        if (descdet->super.myrank == root){
            *determinant = *((double*)(descdet->super.data_of(&descdet->super, 0, 0)->device_copies[0]->device_private));
        }
        MPI_Bcast((void*)determinant, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
        logdet = 2 * (*determinant);
        STOP_TIMING(logdet_calculate);
        VERBOSE(" Done.\n");

#if DEBUG_INFO
        sum = ParsecDZSum(parsec, (parsec_tiled_matrix_dc_t *)descZ);
        if (rank == 0)
            fprintf(stderr, "\nZ: sum before dtrsm %.17g\n", sum);
#endif

        VERBOSE("Solving the linear system ...");
        START_TIMING(time_solve);
        parsec_taskpool_t *parsec_dtrsm = dplasma_dtrsm_New(dplasmaLeft, uplo, dplasmaNoTrans, dplasmaNonUnit,
                                                             1.0, (parsec_tiled_matrix_t*)(&hicma_data->dcA0),
                                                             descZ);
#if NO_GPU_DPLASMA
        disable_GPU(parsec_dtrsm);
#endif
        STOP_TIMING(time_solve);
        flops += FLOPS_DTRSM(PlasmaLeft, N, NRHS);
        VERBOSE(" Done.\n");

#if DEBUG_INFO
        sum = ParsecDZSum(parsec, (parsec_tiled_matrix_dc_t *)descZ);
        if (rank == 0)
            fprintf(stderr, "\nZ: sum after dtrsm %.17g\n", sum);
#endif

#if OWN_DOT
        parsec_dmatrix_dot(parsec, descZ, descZ, descproduct);
#else
        VERBOSE("Calculating the MLE likelihood function ...");
        START_TIMING(time_mle);
        parsec_matrix_block_cyclic_t descZ0;
        parsec_matrix_block_cyclic_init(&descZ0, PARSEC_MATRIX_DOUBLE, PARSEC_MATRIX_TILE,
                                         rank, descZ->mb, descZ->nb, descZ->m,
                                         1, 0, 0, descZ->m, 1, nodes, 1, 1, 1, 0, 0);
        descZ0.mat = calloc((size_t)descZ0.super.nb_local_tiles *
				(size_t)descZ0.super.bsiz,
				(size_t)parsec_datadist_getsizeoftype(descZ0.super.mtype));
		parsec_data_collection_set_key((parsec_data_collection_t*)&descZ0, "descZ0");
        dplasma_dlacpy( parsec, matrix_UpperLower, descZ, (parsec_tiled_matrix_t *)&descZ0 );
		parsec_taskpool_t *parsec_dgemm = NULL;
		parsec_dgemm = dplasma_dgemm_New( PlasmaTrans, PlasmaNoTrans, 1.0, descZ, (parsec_tiled_matrix_t *)&descZ0, 0.0, descproduct );

#if NO_GPU_DPLASMA
        disable_GPU(parsec_dgemm);
#endif
        if ( parsec_dgemm != NULL )
		{
			parsec_context_add_taskpool( parsec, (parsec_taskpool_t*)parsec_dgemm);
			parsec_context_start( parsec );
			parsec_context_wait( parsec );
			dplasma_dgemm_Destruct( parsec_dgemm );
		}

		parsec_data_free( descZ0.mat );
		parsec_tiled_matrix_destroy( (parsec_tiled_matrix_t*)&descZ0 );
#endif

        root = descproduct->super.rank_of(&descproduct->super, 0, 0);
        if (descproduct->super.myrank == root)
            *product = *((double*)(descproduct->super.data_of(&descproduct->super, 0, 0)->device_copies[0]->device_private));
        MPI_Bcast((void*)product, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
        // remove if/else, no comaprison of kernel name with "mattern"
        loglik = -0.5 * (*product) - 0.5 * logdet - (double) (N / 2.0) * log(2.0 * PI);
        variance = theta[0];
        STOP_TIMING(time_mle);
        VERBOSE(" Done.\n");

#if defined(PARSEC_HAVE_MPI)
        MPI_Bcast(&loglik, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast((void*)theta, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD); //TODO: not in cham implementation,does it work in C?
#endif

#if SWITCH_TO_DP
/* parameters to swith back to dp for mix-precision */
int switch_to_dp = 0;
double loglik_pre = 0.0;
        if (fabs(loglik_pre - loglik) < threshold_to_dp)
            switch_to_dp = 1;
        if (rank == 0 && data->precision != 0)
            fprintf(stderr, "SWITCH_TO_DP iter %d : %d -- %lf %lf\n", iter_count, switch_to_dp, fabs(loglik_pre - loglik), threshold_to_dp);
#endif

#if defined(PARSEC_HAVE_MPI)
    if (rank == 0)
#endif
    {
        fprintf(stderr, "\n------ddotproduct: %.17lf ", *product);
        fprintf(stderr, "\n------logdet: %.17lf ", logdet);
        if (num_params == 3) {
            fprintf(stderr, " %3d- Model Parameters (variance, range, smoothness): (%.10lf, %.10lf, %.10lf) ----> LogLi: %.17lf\n",
                    iter_count, theta[0], theta[1], theta[2], loglik);
        } else if (num_params == 7) {
            if (rank == 0)
                fprintf(stderr, "Optimization Start: %3d- Model Parameters (theta0, theta1, theta2, theta3, theta4, theta5, theta6): (%.10lf, %.10lf, %.10lf, %.10lf, %.10lf, %.10lf, %.10lf)  ----> LogLi: %.17lf\n\n",
                        iter_count, theta[0], theta[1], theta[2], theta[3], theta[4], theta[5], theta[6], loglik);
        } else {
            fprintf(stderr, " %3d- Model Parameters (variance, range, smoothness, noise): (%.10lf, %.10lf, %.10lf, %.10lf) ----> LogLi: %.17lf\n",
                    iter_count, theta[0], theta[1], theta[2], theta[3], loglik);
        }
        // if (aConfigurations.GetLogger == 1)
        //     fprintf(aConfigurations.GetLoggerPath, " %3d- Model Parameters (variance, range, smoothness): (%.10lf, %.10lf, %.10lf) ----> LogLi: %.17lf\n",
        //             iter_count, theta[0], theta[1], theta[2], loglik);

        fprintf(stderr, " ---- Facto Time: %.2lf\n", hicma_params->time_hicma);
        fprintf(stderr, " ---- logdet Time: %.2lf\n", logdet_calculate);
        fprintf(stderr, " ---- dtrsm Time: %.2lf\n", time_solve);
        fprintf(stderr, " ---- MLE Time: %.2lf\n", time_mle);
        fprintf(stderr, " ---- Matrix Generation Time: %.2lf\n", matrix_gen_time);
        fprintf(stderr, " ---- Total Time: %.2lf\n", matrix_gen_time + hicma_params->time_hicma + logdet_calculate + time_solve);
        fprintf(stderr, " ---- Gflop/s: %.2lf\n",
                ((time_facto + time_solve) == 0.0) ? 0.0 : flops / 1e9 / (time_facto + time_solve));
        fprintf(stderr, "\n\n");
    }

 // for experiments and benchmarking
    accumulated_executed_time = Results::GetInstance()->GetTotalModelingExecutionTime() + time_facto + logdet_calculate + time_solve;
    Results::GetInstance()->SetTotalModelingExecutionTime(accumulated_executed_time);
    accumulated_flops = Results::GetInstance()->GetTotalModelingFlops() + (flops / 1e9 / (time_facto + time_solve));
    Results::GetInstance()->SetTotalModelingFlops(accumulated_flops);
    Results::GetInstance()->SetMLEIterations(iter_count + 1);
    // Results::GetInstance()->SetMaximumTheta(vector<double>(theta, theta + num_params));
    Results::GetInstance()->SetLogLikValue(loglik);
#if SWITCH_TO_DP
    if (fabs(loglik_pre - loglik) < threshold) {
        if (rank == 0)
            // fprintf(stderr, "\nForce terminate iter %d : %.17lf %.17lf %p\n\n", iter_count - 1,
            //         fabs(loglik_pre - loglik), threshold, *data->opt);
        nlopt_force_stop(*data->opt);
    }
    if (1 || loglik_pre < loglik || loglik_pre == 0.0) {
        if (rank == 0)
            fprintf(stderr, "\nSet loglik_pre %d : %lf %lf\n\n", iter_count - 1, loglik, loglik_pre);
        loglik_pre = loglik;
    }
#endif

    return loglik;
}