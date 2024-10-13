
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

    return 0;
}