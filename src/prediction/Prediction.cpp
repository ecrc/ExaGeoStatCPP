
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file Prediction.cpp
 * @brief Contains the implementation of the Prediction class.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2024-02-04
**/

#include <cstring>

#ifdef USE_MKL
#include <mkl_service.h>
#endif

#include <prediction/Prediction.hpp>
#include <prediction/PredictionHelpers.hpp>
#include <prediction/PredictionAuxiliaryFunctions.hpp>
#include <linear-algebra-solvers/LinearAlgebraFactory.hpp>

#if !DEFAULT_RUNTIME
#include <runtime/parsec/ParsecHeader.h>
extern "C"{
#include <runtime/parsec/jdf/JobDescriptionFormat.h>
}
#endif

using namespace std;

using namespace exageostat::prediction;
using namespace exageostat::dataunits;
using namespace exageostat::results;
using namespace exageostat::configurations;
using namespace exageostat::common;

template<typename T>
void Prediction<T>::PredictMissingData(unique_ptr<ExaGeoStatData<T>> &aData, Configurations &aConfigurations,
                                       T *apMeasurementsMatrix, const kernels::Kernel<T> &aKernel,
                                       Locations<T> *apTrainLocations, Locations<T> *apTestLocations) {                            
#if !DEFAULT_RUNTIME
    double *Zobs;
    double *z_actual;
    double *Zmiss;
    int p_grid = aConfigurations.GetPGrid();
    int q_grid = aConfigurations.GetQGrid();
    int dts = aConfigurations.GetDenseTileSize();
    int mse_flag = 1;
    int i, j, number_of_mspe = 3;
    int p = aKernel.GetVariablesNumber();
    int z_miss_number, n_z_obs;
    bool can_predict = true;
    int num_params = aKernel.GetParametersNumbers();
    int precision = aConfigurations.GetPrecision();
	int nodes = aConfigurations.GetCoresNumber();
	int rank = aConfigurations.GetMaxRank();
    int N = aConfigurations.GetProblemSize();

	parsec_matrix_uplo_t uplo = PARSEC_MATRIX_LOWER;

    for (i = 0; i < num_params; i++) {
        if (aConfigurations.GetEstimatedTheta()[i] == -1) {
            can_predict = false;
            break;
        }
    }

    if (!can_predict &&
        (aConfigurations.GetIsMLOEMMOM() || aConfigurations.GetIsMSPE() || aConfigurations.GetIsFisher())) {
        throw runtime_error(
                "Can't predict without an estimated theta, please either pass --etheta or run the modeling module before prediction");
    }
    
    if (!apTrainLocations && !apTestLocations) {
        z_miss_number = aConfigurations.GetUnknownObservationsNb();
        n_z_obs = aConfigurations.CalculateZObsNumber();
        z_actual = new double[z_miss_number * p];
    } else {
        z_miss_number = apTestLocations->GetSize();
        n_z_obs = apTrainLocations->GetSize();
        z_actual = nullptr;
    }
   
    aConfigurations.SetObservationNumber(n_z_obs);

    VERBOSE("\t- Total number of Z: " << aConfigurations.GetProblemSize())
    LOGGER("\t- Number of Z Miss: " << z_miss_number)
    LOGGER("\t- Number of Z observations: " << n_z_obs)

    //memory allocation
    Zobs 	= (double *) malloc(n_z_obs * sizeof(double));
    Zmiss	= (double *) malloc(z_miss_number * sizeof(double));

    // Initialize prediction descriptors
	if(z_miss_number <= 0)
	{
		fprintf(stderr," Number of missing values should be positive value\n");
		return;
	}

    //bi-variate case
    auto kernel_fun = aConfigurations.GetKernelName();
	if(     kernel_fun == "bivariate_matern_parsimonious" ||
			kernel_fun == "bivariate_matern_parsimonious_profile" ||
			kernel_fun == "bivariate_matern_parsimonious2"||
			kernel_fun == "bivariate_matern_parsimonious2_profile" ||
			kernel_fun == "bivariate_matern_flexible" ||
			kernel_fun == "bivariate_matern_flexible_profile" ||
			kernel_fun == "bivariate_matern_flexible2" ||
			kernel_fun == "bivariate_matern_flexible2_profile" )
	{
		n_z_obs*=2;
		z_miss_number*=2;
	}

    aData->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_Z_MISS);
    aData->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_C11);
    aData->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_C12);
    aData->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_C21);
    aData->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_C22);
    aData->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_Z_TRACE);
    aData->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_Z_OBS);

    aData->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_MSE);
    aData->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_MSE_1);
    aData->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_MSE_2);
    aData->GetDescriptorData()->SetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_Z_Actual);

    parsec_matrix_block_cyclic_t *descZmiss = (aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_Z).parsec_desc);
    parsec_matrix_block_cyclic_t *descC11 = (aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_C11).parsec_desc);
    parsec_matrix_block_cyclic_t* descC21 = (aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_C21).parsec_desc);
    parsec_matrix_block_cyclic_t *descC12 = (aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_C12).parsec_desc);
    parsec_matrix_sym_block_cyclic_t* descC22 = new parsec_matrix_sym_block_cyclic_t;//(aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_C22).parsec_desc);
    parsec_matrix_block_cyclic_t *descZtrace = (aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_Z_TRACE).parsec_desc);
    parsec_matrix_block_cyclic_t *descZobs = (aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_Z_OBS).parsec_desc);

	parsec_matrix_block_cyclic_init(descZobs, PARSEC_MATRIX_DOUBLE, PARSEC_MATRIX_TILE,
			rank, dts, dts, n_z_obs, 1, 0, 0,
			n_z_obs, 1, 1, 1, 1, 1, 0, 0);
	descZobs->mat = calloc((size_t)descZobs->super.nb_local_tiles *
			(size_t)descZobs->super.bsiz,
			(size_t)parsec_datadist_getsizeoftype(descZobs->super.mtype));
	parsec_data_collection_set_key((parsec_data_collection_t*)descZobs, "descZobs");

	parsec_matrix_block_cyclic_init(descZmiss, PARSEC_MATRIX_DOUBLE, PARSEC_MATRIX_TILE,
			rank, dts, dts, z_miss_number, 1, 0, 0,
			z_miss_number, 1, nodes, 1, 1, 1, 0, 0);
	descZmiss->mat = calloc((size_t)descZmiss->super.nb_local_tiles *
			(size_t)descZmiss->super.bsiz,
			(size_t)parsec_datadist_getsizeoftype(descZmiss->super.mtype));
	parsec_data_collection_set_key((parsec_data_collection_t*)descZmiss, "descZmiss");

	parsec_matrix_block_cyclic_init(descC12, PARSEC_MATRIX_DOUBLE, PARSEC_MATRIX_TILE,
			rank, dts, dts, z_miss_number, n_z_obs, 0, 0,
			z_miss_number, n_z_obs, p_grid, nodes/p_grid, 1, 1, 0, 0);
	descC12->mat = calloc((size_t)descC12->super.nb_local_tiles *
			(size_t)descC12->super.bsiz,
			(size_t)parsec_datadist_getsizeoftype(descC12->super.mtype));
	parsec_data_collection_set_key((parsec_data_collection_t*)descC12, "descC12");

	parsec_matrix_sym_block_cyclic_init(descC22, PARSEC_MATRIX_DOUBLE,
			rank, dts, dts,z_miss_number, z_miss_number, 0, 0,
			z_miss_number, z_miss_number, p_grid,nodes/p_grid ,uplo);
	descC22->mat = calloc((size_t)descC22->super.nb_local_tiles *
			(size_t)descC22->super.bsiz,
			(size_t)parsec_datadist_getsizeoftype(descC22->super.mtype));
	parsec_data_collection_set_key((parsec_data_collection_t*)descC22, "descC22");

	parsec_matrix_block_cyclic_init(descZtrace, PARSEC_MATRIX_DOUBLE, PARSEC_MATRIX_TILE,
			rank, dts, dts, z_miss_number, 1, 0, 0,
			z_miss_number, 1, p_grid, 1, 1, 1, 0, 0);
	descZtrace->mat = calloc((size_t)descZtrace->super.nb_local_tiles *
			(size_t)descZtrace->super.bsiz,
			(size_t)parsec_datadist_getsizeoftype(descZtrace->super.mtype));
	parsec_data_collection_set_key((parsec_data_collection_t*)descZtrace, "descZtrace");

	parsec_matrix_block_cyclic_init(descC11, PARSEC_MATRIX_DOUBLE, PARSEC_MATRIX_TILE,
			rank, dts, dts, z_miss_number, z_miss_number, 0, 0,
			z_miss_number, z_miss_number, p_grid, nodes/p_grid, 1, 1, 0, 0);
	descC11->mat = calloc((size_t)descC11->super.nb_local_tiles *
			(size_t)descC11->super.bsiz,
			(size_t)parsec_datadist_getsizeoftype(descC11->super.mtype));
	parsec_data_collection_set_key((parsec_data_collection_t*)descC11, "descC11");

	parsec_matrix_block_cyclic_init(descC21, PARSEC_MATRIX_DOUBLE, PARSEC_MATRIX_TILE,
			rank, dts, dts, n_z_obs, z_miss_number, 0, 0,
			n_z_obs, z_miss_number, p_grid, nodes/p_grid, 1, 1, 0, 0);
	descC21->mat = calloc((size_t)descC21->super.nb_local_tiles *
			(size_t)descC21->super.bsiz,
			(size_t)parsec_datadist_getsizeoftype(descC21->super.mtype));
	parsec_data_collection_set_key((parsec_data_collection_t*)descC21, "descC21");

	parsec_matrix_block_cyclic_t *descmse = (aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_MSE).parsec_desc);
	parsec_matrix_block_cyclic_t *descmse1 = (aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_MSE_1).parsec_desc);
	parsec_matrix_block_cyclic_t *descmse2 = (aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_MSE_2).parsec_desc);
	parsec_matrix_block_cyclic_t *descZactual = (aData->GetDescriptorData()->GetDescriptor(PARSEC_DESCRIPTOR, DESCRIPTOR_Z_Actual).parsec_desc);

	if( mse_flag == 1) {
		parsec_matrix_block_cyclic_init(descZactual, PARSEC_MATRIX_DOUBLE, PARSEC_MATRIX_TILE,
				rank, dts, dts, z_miss_number, 1, 0, 0,
				z_miss_number, 1, nodes, 1, 1, 1, 0, 0);
		descZactual->mat = calloc((size_t)descZactual->super.nb_local_tiles *
				(size_t)descZactual->super.bsiz,
				(size_t)parsec_datadist_getsizeoftype(descZactual->super.mtype));
		parsec_data_collection_set_key((parsec_data_collection_t*)descZactual, "descZactual");

		parsec_matrix_block_cyclic_init(descmse, PARSEC_MATRIX_DOUBLE, PARSEC_MATRIX_TILE,
				rank, dts, dts, 1, 1, 0, 0,
				1, 1, nodes, 1, 1, 1, 0, 0);
		descmse->mat = calloc((size_t)descmse->super.nb_local_tiles *
				(size_t)descmse->super.bsiz,
				(size_t)parsec_datadist_getsizeoftype(descmse->super.mtype));
		parsec_data_collection_set_key((parsec_data_collection_t*)descmse, "descmse");
	}

    // temporary fetch locations from exageostat data object
    auto locations = aData->GetLocations();
    location l1;
    l1.x = (double *) locations->GetLocationX();
    l1.y = (double *) locations->GetLocationY();
    l1.z = (double *) locations->GetLocationZ();    
    double avg_pred_value=0.0;
    int pred_samples = 1;
    j = 0;
    auto *parsec = (parsec_context_t *) ExaGeoStatHardware::GetParsecContext();

    for (j = 0; j < pred_samples; j++)
    {
        // pick random points 
        /* Initialization */
            location l;
            location* lmiss = new location();
            location* lobs= new location();
            location* lm = new location();
            double *Z;                      
            int i = 0;
        
            /* Memory allocation */
            Z       = (double *) malloc(N * sizeof(double));
            l.x     = (double *) malloc(N * sizeof(double));
            l.y     = (double *) malloc(N * sizeof(double)); 
        
            if(l1.z != NULL)
                l.z     = (double *) malloc(N * sizeof(double));
            else
                l.z = NULL;
        
            parsec_tiled_matrix_t *descZcpy = &((aData->GetDescriptorData()->GetDescriptor(DescriptorType::PARSEC_DESCRIPTOR, DescriptorName::DESCRIPTOR_Z_COPY).parsec_desc)->super);
            assert( descZcpy->m == N );
                
            /* Copy observed measurments */
            ParsecGetZObs( parsec, descZcpy, Z, N );
        
            /* Broadcast Z */
            int root = descZcpy->super.rank_of(&descZcpy->super, 0, 0);
            MPI_Bcast( Z, N, MPI_DOUBLE, root, MPI_COMM_WORLD );
        
            for( i = 0; i < N ; i++)        
            {
                l.x[i]=l1.x[i];
                l.y[i]=l1.y[i];
        
                if(l1.z != NULL)
                    l.z[i]=l1.z[i];
            }               
        
            // print_dmatrix("l.x", N, 1, l.x, N);         
            // print_dmatrix("l.y", N, 1, l.y, N);

            printf("N %d z_miss_number %d n_z_obs %d\n", N, z_miss_number, n_z_obs);
            // shuffle(Z, &l, N);-> replaced with its body

            if (N > 1) {
                printf("inside shuffle\n");
                size_t i;
                for (i = 0; i < N - 1; i++) 
                {
                    size_t j = i + rand() / (RAND_MAX / (N - i) + 1);

                    double t = Z[j];
                    Z[j] = Z[i];
                    Z[i] = t;
                    double xtemp = l.x[j];
                    l.x[j] = l.x[i];
                    l.x[i] = xtemp;
                    double ytemp = l.y[j];
                    l.y[j] = l.y[i];
                    l.y[i] = ytemp;
                }
            }

            //print_dmatrix("Z after", N, 1, Z, N);
        
            for( i = 0; i < n_z_obs ; i++)
                Zobs[i] = Z[z_miss_number+i];
        
            for ( i = 0; i < z_miss_number; i++)
                z_actual[i] = Z[i];

            //print_datrix("Zobs", n_z_obs, 1, Zobs, n_z_obs);
            //print_dmatrix("z_actual", z_miss_number, 1, z_actual, z_miss_number);
        
            // &lmiss = &(data->&lmiss);
            // lobs = &(data->lobs); ARE THEY SET DIFFERENTLY FEL CHAM??

            lmiss->x = l.x;
            lmiss->y = l.y;

            lobs->x  = &l.x[z_miss_number];
            lobs->y  = &l.y[z_miss_number];

            if(l1.z != NULL)
            {
                lmiss->z = l.z;
                lobs->z  = &l.z[z_miss_number];
            }

            //should be replaced with already existing sorting
            PredictionHelpers<T>::LocationsZSortInPlace(n_z_obs, lobs, Zobs);
            PredictionHelpers<T>::LocationsZSortInPlace(z_miss_number, lmiss, z_actual);
           
            double time_solve = 0.0;
            double mat_gen_time = 0.0;
            double time_gemm = 0.0;
            double time_mse = 0.0;
            double flops = 0.0;
            int info = 0;
       
            parsec_tiled_matrix_t *descdet = &((aData->GetDescriptorData()->GetDescriptor(DescriptorType::PARSEC_DESCRIPTOR, DescriptorName::DESCRIPTOR_DETERMINANT).parsec_desc)->super);
            int mserror                   = 0.0;
        
            if(z_miss_number <= 0)
            {
                fprintf(stderr," Number of missing values should be positive value\n");
                return ;
            }
        
            int HNB = aConfigurations.GetHNB();
            int success;
            double *variance;
            /* Synchronization before prediction */
            MPI_Barrier( MPI_COMM_WORLD );
        
            VERBOSE("Copy measurments vector to descZobs descriptor...");
            LapackToTile( parsec, &(descZobs->super), Zobs, n_z_obs );
            VERBOSE(" Done.\n");
        
            if( z_actual != NULL)
            {
                printf("parsec context: %p, descZactual: %p, Zactual: %p, nZmiss: %d\n",parsec, &(descZactual->super), z_actual, z_miss_number);
                VERBOSE("Copy actual measurments vector to descZactual descriptor...");
                LapackToTile( parsec, &(descZactual->super), z_actual, z_miss_number );
                VERBOSE(" Done.\n");
            }
        
        // #if defined(PARSEC_HAVE_MPI)
        //     MPI_Bcast(variance, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
        // #endif
        auto initial_theta = aConfigurations.GetInitialTheta();
        int distance_matrix = aConfigurations.GetDistanceMetric();

        
            printf("estimated parameters:");
            i = 0;
            if( 0 == rank )
            {
                for(i=0; i<num_params; i++)
                {
                    printf("%.8f,", initial_theta[i]);
        
                }
                printf(")\n");
            }
                
        #if defined(EXAGEOSTAT_USE_CUDA) 
            /* Reset cache on GPU before cholesky */
            parsec_devices_release_memory();
            parsec_devices_reset_load(parsec);
        #endif
        
            /* Generate C22 covariance matrix */
            START_TIMING(mat_gen_time);
            VERBOSE("Generate C22 Covariance Matrix... (Prediction Stage)");
            ParsecDMatrixGeneration( parsec, &(descC22->super), lobs, lobs,
                    lm, initial_theta.data(), distance_matrix,"univariate_matern_stationary", descC22->super.lmt );
            VERBOSE(" Done.\n");
        
            /* Generate C12 covariance matrix */
            VERBOSE("Generate C12 Covariance Matrix... (Prediction Stage)");
            ParsecDMatrixGeneration( parsec, &(descC12->super), lmiss, lobs,
            lm, initial_theta.data(), distance_matrix,"univariate_matern_stationary", descC12->super.lmt );
            VERBOSE(" Done.\n");
        
            /* Generate C11 covariance matrix */
            VERBOSE("Generate C11 Covariance Matrix... (Prediction Stage)");
            ParsecDMatrixGeneration( parsec, &(descC12->super), lmiss, lmiss,
            lm, initial_theta.data(), distance_matrix,"univariate_matern_stationary", descC12->super.lmt );
            VERBOSE(" Done.\n");
        
            /* Generate C21 covariance matrix */
            VERBOSE("Generate C21 Covariance Matrix... (Prediction Stage)");
            ParsecDMatrixGeneration( parsec, &(descC21->super), lobs, lmiss,
            lm, initial_theta.data(), distance_matrix,"univariate_matern_stationary", descC21->super.lmt );
            VERBOSE(" Done.\n");
        
            STOP_TIMING(mat_gen_time);
        
            /* Start prediction */
            VERBOSE("Calculate dposv C22 Covariance Matrix... (Prediction Stage)");
                
            /* Define band_size */
            if( 0 == aConfigurations.GetPrecision() ) {
                int NT = descC22->super.lmt ;
                aConfigurations.SetDenseBandDP(NT); 
                aConfigurations.SetBandDenseSP(NT);
            }
        
            VERBOSE("Calculate hsdpotrf C22 Covariance Matrix... (Prediction Stage)");
            /* Start prediction */
            START_TIMING(time_solve);
            success = dplasma_dpotrf(parsec, uplo, &(descC22->super));
            VERBOSE(" Done.\n");
        
            if( 0 != success && 0 == rank ) {
                fprintf(stderr, "\n\nFactorization cannot be performed..\n The matrix is not positive definite: %d\n", success);
                exit(1);
            }
        
            VERBOSE("Triangular solve 1... (Prediction Stage)");
            dplasma_dtrsm( parsec, dplasmaLeft, dplasmaLower, dplasmaNoTrans, dplasmaNonUnit, 1.0, &(descC22->super), &(descZobs->super) );
            VERBOSE(" Done.\n");
        
            VERBOSE("Triangular solve 2... (Prediction Stage)");
            dplasma_dtrsm( parsec, dplasmaLeft, dplasmaLower, dplasmaTrans, dplasmaNonUnit, 1.0, &(descC22->super), &(descZobs->super) );
            VERBOSE(" Done.\n");
        
            STOP_TIMING(time_solve);
        
            flops = flops + FLOPS_DPOTRF(n_z_obs);
            flops = flops + 2 * FLOPS_DTRSM(dplasmaLeft, n_z_obs, n_z_obs);
        
            START_TIMING(time_gemm);
            VERBOSE("Calculate dgemm Zmiss= C12 * Zobs Covariance Matrix... (Prediction Stage)");   
            parsec_taskpool_t *parsec_dgemm = NULL;
            parsec_dgemm = dplasma_dgemm_New( dplasmaNoTrans, dplasmaNoTrans, 1.0, &(descC12->super), &(descZobs->super), 0.0, &(descZmiss->super) ); 
        #if NO_GPU_DPLASMA 
            // disable_GPU( parsec_dgemm );
            #if defined(EXAGEOSTAT_USE_CUDA)
            for (int i = 0; i < parsec_nb_devices; i++) {
                parsec_device_module_t *device = parsec_mca_device_get(i);
                if( device->type == PARSEC_DEV_CUDA )
                parsec_dgemm->devices_index_mask &= ~(1<<i);
            }
            #endif
        #endif
            if ( parsec_dgemm != NULL )
            {
                parsec_context_add_taskpool( parsec, (parsec_taskpool_t*)parsec_dgemm);
                parsec_context_start( parsec );
                parsec_context_wait( parsec );
                dplasma_dgemm_Destruct( parsec_dgemm );
            }
            flops = flops + FLOPS_DGEMM(z_miss_number, n_z_obs, n_z_obs);
            VERBOSE(" Done.\n");
            STOP_TIMING(time_gemm);
        
            /* Schur_complement */
            VERBOSE("Calculate dposv C22 Covariance Matrix and C21...(Prediction Stage - Schur_complement)");
            parsec_taskpool_t * parsec_dtrsm= NULL;
            parsec_dtrsm = dplasma_dtrsm_New( dplasmaLeft, dplasmaLower, dplasmaNoTrans, dplasmaNonUnit, 1.0, &(descC22->super), &(descC21->super) );
        #if NO_GPU_DPLASMA 
            disable_GPU( parsec_dtrsm );
        #endif
            if ( parsec_dtrsm != NULL )
            {
                parsec_context_add_taskpool( parsec, parsec_dtrsm );
                parsec_context_start( parsec );
                parsec_context_wait( parsec );
                dplasma_dtrsm_Destruct( parsec_dtrsm );
            }
            VERBOSE(" Done.\n");
        
            VERBOSE("Calculate dposv C22 Covariance Matrix and C21...(Prediction Stage - Schur_complement)");
            parsec_dtrsm = dplasma_dtrsm_New( dplasmaLeft, dplasmaLower, dplasmaTrans, dplasmaNonUnit, 1.0, &(descC22->super), &(descC21->super) );
        #if NO_GPU_DPLASMA 
            disable_GPU( parsec_dtrsm );
        #endif
            if ( parsec_dtrsm != NULL )
            {
                parsec_context_add_taskpool( parsec, parsec_dtrsm );
                parsec_context_start( parsec );
                parsec_context_wait( parsec );
                dplasma_dtrsm_Destruct( parsec_dtrsm );
            }
        
            flops = flops + FLOPS_DTRSM(dplasmaLeft, n_z_obs, n_z_obs);
            VERBOSE(" Done.\n");
        
            VERBOSE("Calculate dgemm C11 = C12 * C21 Covariance Matrix... (Prediction Stage - Schur_complement)");
            dplasma_dgemm( parsec, dplasmaNoTrans, dplasmaNoTrans, -1.0, &(descC12->super), &(descC21->super), 1.0, &(descC11->super) );
            flops = flops + FLOPS_DGEMM(z_miss_number, n_z_obs, n_z_obs);
            VERBOSE(" Done.\n");
        
            /* Use data->det to store the trace summation */
            int det=0;
            T *determinant = aData->GetDescriptorData()->GetDescriptorMatrix(PARSEC_DESCRIPTOR, DESCRIPTOR_DETERMINANT);
            *determinant = 0;

            VERBOSE("Calculate trace estimation... (Prediction Stage - Schur_complement)");
            MLEDtrace(parsec, &(descC11->super), (descdet), &(descZtrace->super));

            root = descdet->super.rank_of(&descdet->super, 0, 0);
            if( descdet->super.myrank == root )
                *determinant = *((double*)(descdet->super.data_of(&descdet->super, 0, 0)->device_copies[0]->device_private));
        #if defined(PARSEC_HAVE_MPI)
            MPI_Bcast( (void *)determinant, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
        #endif
            VERBOSE(" Done.\n");
        
            double mean = (*determinant)/z_miss_number;
            double sd = 0;
            double sum_trace = 0.0;
            double* Ztrace = (double *) malloc(z_miss_number  * sizeof(double));
            ParsecGetZObs( parsec, &(descZtrace->super), Ztrace, z_miss_number ); 
        
        #if defined(PARSEC_HAVE_MPI)
            MPI_Bcast( (void *)Ztrace, z_miss_number, MPI_DOUBLE, root, MPI_COMM_WORLD);
        #endif
        
            for(int i = 0; i < z_miss_number; i++) {
                sd += pow((Ztrace[i] - mean), 2);
                sum_trace += Ztrace[i];
            }
            sd = sqrt(sd/z_miss_number);
            free( Ztrace );
        
            // data->trace_sum = sum_trace;
            // data->trace_mean = mean;
            // data->trace_sd = sd;
        
        #if defined(PARSEC_HAVE_MPI)
            if(descC22->super.super.myrank == 0)
        #endif
            {
                fprintf(stderr, "Trace_sum %d : %le %le\n", rank, *determinant, sum_trace);
                printf("Trace estimation (trace standard deviation ): %6.4e\n", sd);
                printf("Trace estimation (trace mean ): %6.4e\n", mean);
                printf("Trace estimation (trace sum ): %6.4e\n", *determinant);
            }
        
            /*****************************************************************************************************/
        
            /* return back descZmiss to zmiss vector */
            ParsecGetZObs( parsec, &(descZmiss->super), Zmiss, z_miss_number ); 
        
            root = descZmiss->super.super.rank_of(&(descZmiss->super.super), 0, 0);
            MPI_Bcast( Zmiss, z_miss_number, MPI_DOUBLE, root, MPI_COMM_WORLD );
        
            /* Estimate Mean Square Error */
            if( z_actual != NULL)
            {
                START_TIMING(time_mse);
                VERBOSE("Calculate Mean Square Error (MSE) ... (Prediction Stage) \n");
                MSECalculation(parsec, &(descZactual->super), &(descZmiss->super), &(descmse->super));
                /* Broadcast determinant */
                root = descmse->super.super.rank_of(&descmse->super.super, 0, 0);
                if( descmse->super.super.myrank == root )
                    mserror = *( (double *)((descmse->super.super.data_of(&descmse->super.super, 0, 0))->device_copies[0]->device_private) );
        
        #if defined(PARSEC_HAVE_MPI)
                MPI_Bcast( (void *)&mserror, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
        #endif
                VERBOSE(" Done.\n");    
                STOP_TIMING(time_mse);
                mserror /= z_miss_number;
            }
            else
                mserror = -1;
        
        if( 0 == rank ) {
            fprintf(stderr,"Prediction Error: %lf \n", mserror);
        }

        avg_pred_value +=mserror;
    }

    //free memory
    free(Zobs);
    free(z_actual);
    free(Zmiss);
 #else
    int i, j;
    bool can_predict = true;
    int num_params = aKernel.GetParametersNumbers();
    for (i = 0; i < num_params; i++) {
        if (aConfigurations.GetEstimatedTheta()[i] == -1) {
            can_predict = false;
            break;
        }
    }

    if (!can_predict &&
        (aConfigurations.GetIsMLOEMMOM() || aConfigurations.GetIsMSPE() || aConfigurations.GetIsFisher())) {
        throw runtime_error(
                "Can't predict without an estimated theta, please either pass --etheta or run the modeling module before prediction");
    }

    int number_of_mspe = 3;
    int p = aKernel.GetVariablesNumber();
    int z_miss_number, n_z_obs;
    T *z_actual;
    if (!apTrainLocations && !apTestLocations) {
        z_miss_number = aConfigurations.GetUnknownObservationsNb();
        n_z_obs = aConfigurations.CalculateZObsNumber();
        z_actual = new T[z_miss_number * p];
    } else {
        z_miss_number = apTestLocations->GetSize();
        n_z_obs = apTrainLocations->GetSize();
        z_actual = nullptr;
    }
    aConfigurations.SetObservationNumber(n_z_obs);
    auto linear_algebra_solver = linearAlgebra::LinearAlgebraFactory<T>::CreateLinearAlgebraSolver(common::EXACT_DENSE);

    VERBOSE("\t- Total number of Z: " << aConfigurations.GetProblemSize())
    LOGGER("\t- Number of Z Miss: " << z_miss_number)
    LOGGER("\t- Number of Z observations: " << n_z_obs)

    // FISHER Prediction Function Call
    if (aConfigurations.GetIsFisher()) {
        LOGGER("\t---- Using Prediction Function Fisher ----")
        T *fisher_results;
        fisher_results = linear_algebra_solver->ExaGeoStatFisherTile(aConfigurations, aData,
                                                                     (T *) aConfigurations.GetEstimatedTheta().data(),
                                                                     aKernel);
        vector<double> fisher_vector;
        fisher_vector.reserve(num_params * num_params); // Reserve memory in advance for efficiency
        for (size_t idx = 0; idx < num_params * num_params; ++idx) {
            fisher_vector.push_back(fisher_results[idx]);
        }
        Results::GetInstance()->SetFisherMatrix(fisher_vector);
        VERBOSE("\t\t- Sd of sigma2, alpha, nu: " << sqrt(fisher_results[0]) << " " << sqrt(fisher_results[4]) << " "
                                                  << sqrt(fisher_results[8]))
        VERBOSE("\t\t- CI for sigma2: " << aConfigurations.GetEstimatedTheta()[0] - Q_NORM * sqrt(fisher_results[0])
                                        << " "
                                        << aConfigurations.GetEstimatedTheta()[0] + Q_NORM * sqrt(fisher_results[0]))

        VERBOSE("\t\t- CI for alpha: " << aConfigurations.GetEstimatedTheta()[1] - Q_NORM * sqrt(fisher_results[4])
                                       << " "
                                       << aConfigurations.GetEstimatedTheta()[1] + Q_NORM * sqrt(fisher_results[4]))
        VERBOSE("\t\t- CI for nu: " << aConfigurations.GetEstimatedTheta()[2] - Q_NORM * sqrt(fisher_results[8]) << " "
                                    << aConfigurations.GetEstimatedTheta()[2] + Q_NORM * sqrt(fisher_results[8]))

        VERBOSE("\t\t- Fisher Matrix:")
        for (i = 0; i < num_params; i++) {
            VERBOSE("\t\t ", true)
            for (j = 0; j < num_params; j++) {
                VERBOSE(fisher_results[i * num_params + j],true)
                if (j != num_params - 1) {
                    VERBOSE(", ",true)
                }
            }
            VERBOSE("")
        }
        delete[] fisher_results;
    }

    if (z_miss_number <= 0) {
        return;
    }

    Results::GetInstance()->SetZMiss(z_miss_number);
    aConfigurations.SetUnknownObservationsNb(z_miss_number);
    T *z_obs = new T[n_z_obs * p];
    T *z_miss = new T[z_miss_number];

    vector<T> avg_pred_value(number_of_mspe);
    auto miss_locations = new Locations<T>(z_miss_number, aConfigurations.GetDimension());
    // Prediction is only supported with 2D.
    auto obs_locations = new Locations<T>(n_z_obs, aConfigurations.GetDimension());

    // We Predict date with only Exact computation. This is a pre-request.
    InitializePredictionArguments(aConfigurations, aData, linear_algebra_solver, z_obs, z_actual, *miss_locations,
                                  *obs_locations, apMeasurementsMatrix, p, apTrainLocations, apTestLocations);

    // MLOE MMOM Auxiliary Function Call
    if (aConfigurations.GetIsMLOEMMOM()) {
        LOGGER("---- Using Auxiliary Function MLOE MMOM ----")
        linear_algebra_solver->ExaGeoStatMLETileMLOEMMOM(aConfigurations, aData,
                                                         (T *) aConfigurations.GetInitialTheta().data(),
                                                         (T *) aConfigurations.GetEstimatedTheta().data(),
                                                         *miss_locations, *obs_locations, aKernel);
    }

    // IDW Auxiliary Function Call
    if (aConfigurations.GetIsIDW()) {
        LOGGER("\t---- Using Auxiliary Function IDW ----")
        T *mspe = new T[number_of_mspe];

        if (!z_actual) {
            z_actual = new T[z_miss_number * p];
            memcpy(z_actual, apMeasurementsMatrix + n_z_obs, z_miss_number * sizeof(T));
        }

        PredictionAuxiliaryFunctions<T>::PredictIDW(z_miss, z_actual, z_obs, z_miss_number, n_z_obs, *miss_locations,
                                                    *obs_locations, mspe);

        vector<double> idw_error;
        idw_error.reserve(number_of_mspe); // Reserve memory in advance for efficiency
        for (size_t idx = 0; idx < number_of_mspe; ++idx) {
            idw_error.push_back(mspe[idx]);
        }

        Results::GetInstance()->SetIDWError(idw_error);
        delete[] mspe;
    }

    // MSPE Prediction Function Call
    if (aConfigurations.GetIsMSPE()) {
        LOGGER("\t---- Using Prediction Function MSPE ----")
        T *prediction_error_mspe;
        if (aConfigurations.GetIsNonGaussian()) {
            prediction_error_mspe = linear_algebra_solver->ExaGeoStatMLENonGaussianPredictTile(aData,
                                                                                               (T *) aConfigurations.GetEstimatedTheta().data(),
                                                                                               z_miss_number, n_z_obs,
                                                                                               z_obs, z_actual, z_miss,
                                                                                               aConfigurations,
                                                                                               *miss_locations,
                                                                                               *obs_locations, aKernel);
        } else {
            prediction_error_mspe = linear_algebra_solver->ExaGeoStatMLEPredictTile(aData,
                                                                                    (T *) aConfigurations.GetEstimatedTheta().data(),
                                                                                    z_miss_number, n_z_obs, z_obs,
                                                                                    z_actual, z_miss, aConfigurations,
                                                                                    *miss_locations, *obs_locations,
                                                                                    aKernel);
        }
        for (i = 0; i < number_of_mspe; i++) {
            avg_pred_value[i] += prediction_error_mspe[i];
        }
        vector<double> z_miss_vector;
        z_miss_vector.reserve(z_miss_number); // Reserve memory in advance for efficiency
        for (size_t idx = 0; idx < z_miss_number; ++idx) {
            z_miss_vector.push_back(z_miss[idx]);
        }
        Results::GetInstance()->SetPredictedMissedValues(z_miss_vector);
        if (z_actual) {
            LOGGER("\t\t- Prediction value: " << avg_pred_value[0])
        }
        delete[] prediction_error_mspe;
    }


    // Due to a leak in Chameleon, exactly trsm We had to free the buffer manually.
#ifdef USE_MKL
    mkl_free_buffers();
#endif

    delete[] z_obs;
    delete[] z_miss;
    delete[] z_actual;
    delete miss_locations;
    delete obs_locations;
#endif
}

template<typename T>
void Prediction<T>::InitializePredictionArguments(Configurations &aConfigurations, unique_ptr<ExaGeoStatData<T>> &aData,
                                                  unique_ptr<linearAlgebra::LinearAlgebraMethods<T>> &aLinearAlgebraSolver,
                                                  T *apZObs, T *apZActual, Locations<T> &aMissLocation,
                                                  Locations<T> &aObsLocation, T *apMeasurementsMatrix, const int &aP,
                                                  Locations<T> *apTrainLocations, Locations<T> *apTestLocations) {
#if DEFAULT_RUNTIME
    int full_problem_size = aConfigurations.GetProblemSize() * aP;
    T *z = new T[full_problem_size];

    aLinearAlgebraSolver->ExaGeoStatGetZObs(aConfigurations, z, full_problem_size, *aData->GetDescriptorData(),
                                            apMeasurementsMatrix, aP);

    if (!apTrainLocations && !apTestLocations) {
        PredictionHelpers<T>::PickRandomPoints(aConfigurations, aData, apZObs, apZActual, z, aMissLocation,
                                               aObsLocation, aP);
    } else {
        for (int i = 0; i < apTrainLocations->GetSize(); ++i) {
            aObsLocation.GetLocationX()[i] = apTrainLocations->GetLocationX()[i];
            aObsLocation.GetLocationY()[i] = apTrainLocations->GetLocationY()[i];
        }
        for (int i = 0; i < apTestLocations->GetSize(); ++i) {
            aMissLocation.GetLocationX()[i] = apTestLocations->GetLocationX()[i];
            aMissLocation.GetLocationY()[i] = apTestLocations->GetLocationY()[i];
        }
        memcpy(apZObs, apMeasurementsMatrix, aObsLocation.GetSize() * sizeof(T));
    }
    delete[] z;
#endif
}