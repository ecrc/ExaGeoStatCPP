
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ClimateEmulator.cpp
 * @brief
 * @details
 * @version 2.0.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2024-09-23
**/

#include <hicma_parsec.h>

int main(int argc, char ** argv)
{
    int iparam[IPARAM_SIZEOF] = {0};
    double dparam[DPARAM_SIZEOF];
    char *cparam[CPARAM_SIZEOF];
    hicma_parsec_params_t params;
    starsh_params_t params_kernel;
    hicma_parsec_data_t data;
    hicma_parsec_matrix_analysis_t analysis;

    /* Init */
    parsec_context_t* parsec = hicma_parsec_init( argc, argv, iparam, dparam, cparam, &params, &params_kernel, &data );

    SYNC_TIME_START();
    //gb24_init_serial(&gb);
    SYNC_TIME_PRINT(params.rank, ("gb24_init\n"));

    // Forward SHT
    SYNC_TIME_START();

    // Forward SHT reshape
    SYNC_TIME_START();
    SYNC_TIME_PRINT(params.rank, ("geqsht_forward_reshape\n"));

    // Make sure it's dense and only off_band maters
    if(params.band_size_dist != 0) {
        if( 0 == params.rank ) {
            fprintf(stderr, RED "Fatal error: band_size_dist= %d needs to be 0\n" RESET, params.band_size_dist);
        }
        params.band_size_dist = 0;
        //return 1;
    }

    if(params.band_size_dense < ceil((double)params.N/params.NB)) {
        if( 0 == params.rank ) {
            fprintf(stderr, RED "Fatal error: the matrix needs to be all DENSE\n" RESET);
        }
        return 1;
    }

    /* Generate matrix */
    SYNC_TIME_START();
    SYNC_TIME_PRINT(params.rank, ("Matrix genneration Matrix norm: norm_global= %le\n", params.norm_global));

    // SYRK
    SYNC_TIME_START();
    dplasma_dsyrk(parsec, dplasmaLower, dplasmaNoTrans,
             1.0, (parsec_tiled_matrix_t *)&data.dcA,
             0.0,  (parsec_tiled_matrix_t *)&data.dcA);
    SYNC_TIME_PRINT(params.rank, ("SYRK\n"));

    // Calculate norm
    SYNC_TIME_START();
    SYNC_TIME_PRINT(params.rank, ("Matrix norm: norm_global= %le\n", params.norm_global));

    /* Analyze matrix before Cholesky */
    hicma_parsec_matrix_pre_analysis( parsec, &data, &params, &params_kernel, &analysis );

    /* HiCMA Cholesky */
    for( int i= 0; i < params.nruns; i++ ) {
        hicma_parsec_potrf( parsec, &data, &params, &analysis );
    }

    if( 0 == params.rank && params.info != 0 && 1 == params.nruns ) {
        fprintf(stderr, "-- Factorization is suspicious (info = %d) ! \n", params.info);
        //exit(params.info);
    }

    /* Analyze matrix after Cholesky */
    hicma_parsec_matrix_post_analysis( parsec, &data, &params, &params_kernel, &analysis );

    SYNC_TIME_PRINT(params.rank, ("mse\n"));

    /* Finalize */
    hicma_parsec_fini( parsec, argc, argv, iparam, dparam, cparam, &params, &params_kernel, &data, &analysis );

    return 0;
}
