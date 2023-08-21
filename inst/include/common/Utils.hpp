
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file Utils.hpp
 * @brief This file contains common functions used in ExaGeoStat software package.
 * @details These functions include so far the VERBOSE macro that prints a message
 * to the console if the verbosity setting is set to "verbose mode.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-03-21
**/

#ifndef EXAGEOSTATCPP_UTILS_HPP
#define EXAGEOSTATCPP_UTILS_HPP

#include <iostream>
#include <string>
#include <chrono>
#include <sys/time.h>

#include <common/Definitions.hpp>
#include <configurations/Configurations.hpp>


/**
 * Verbose macro for logging and debugging mode
 */
#define VERBOSE(msg) { if(exageostat::configurations::Configurations::GetRunMode() == exageostat::common::RunMode::VERBOSE_MODE) std::cout << msg << std::endl; }

/**
 * Verbose macro for logging any failure operation.
 */
#define FAILURE_LOGGER(failed, msg) \
    if (failed){ \
        throw std::runtime_error(msg);\
    }

/**
 * FLOPS MACRO that compute the number of floating point multiplications required
 * for the Cholesky factorization of an n-by-n matrix.
 */
#define FMULS_POTRF(__n) ((double)(__n) * (((1. / 6.) * (double)(__n) + 0.5) * (double)(__n) + (1. / 3.)))

/**
 * FLOPS MACRO that compute the number of floating point additions required
 * for the Cholesky factorization of an n-by-n matrix.
 */
#define FADDS_POTRF(__n) ((double)(__n) * (((1. / 6.) * (double)(__n)      ) * (double)(__n) - (1. / 6.)))

/**
 * FLOPS MACRO that combine the above two macros to compute the
 * total number of floating point operations required for the Cholesky factorization of an n-by-n matrix.
 */
#define FLOPS_DPOTRF(__n) (     FMULS_POTRF((__n)) +       FADDS_POTRF((__n)) )

/**
 * FLOPS MACRO that compute the number of floating point multiplications
 * required for a triangular matrix-matrix multiplication operation where one matrix is 2-by-n and the other is n-by-n.
 */
#define FMULS_TRMM_2(__m, __n) (0.5 * (double)(__n) * (double)(__m) * ((double)(__m)+1.))

/**
 * FLOPS MACRO that compute the number of floating point multiplications required
 * for a triangular matrix-matrix multiplication operation where one matrix is m-by-n and the other is n-by-n.
 */
#define FMULS_TRMM(__side, __m, __n) ( ( (__side) == ChamLeft ) ? FMULS_TRMM_2((__m), (__n)) : FMULS_TRMM_2((__n), (__m)) )

/**
 * Aliases for the FMULS_TRMM macro,
 */
#define FMULS_TRSM FMULS_TRMM
#define FADDS_TRSM FMULS_TRMM

/**
 * FLOPS macrp that combine the FMULS_TRMM and FADDS_TRSM macros to compute the total number of floating point operations
 * required for a triangular matrix solve operation where one matrix is m-by-n and the other is n-by-n.
 */
#define FLOPS_DTRSM(__side, __m, __n) (     FMULS_TRSM(__side, (__m), (__n)) +       FADDS_TRSM(__side, (__m), (__n)) )

/**
 * Timing macro to start timing.
 */
#define START_TIMING(t) auto t##_start = std::chrono::high_resolution_clock::now()

/**
 * Timing macro to stop timing.
 */
#define STOP_TIMING(t) auto t##_end = std::chrono::high_resolution_clock::now(); \
                    t = std::chrono::duration_cast<std::chrono::duration<double>>(t##_end - t##_start).count()

/// static bool to make sure that print summary is only performed once.
static bool is_printed = false;

/**
 * @brief print the summary of MLE inputs.
 * @param[in] N Number of Locations
 * @param[in] ncores Number of Threads per node
 * @param[in] gpus GPUs
 * @param[in] ts Dense Tile Size
 * @param[in] computation
 * @param[in] p_grid
 * @param[in] q_grid
 * @param[in] precision Double or Single Precision
 * @return void
 */
inline void
PrintSummary(int N, int ncores, int gpus, int ts, int p_grid, int q_grid,
             exageostat::common::Precision precision) {
    if (!is_printed) {
#if defined(CHAMELEON_USE_MPI)
        if ( MORSE_My_Mpi_Rank() == 0 )
        {
#endif
        fprintf(stderr, "********************SUMMARY**********************\n");
        fprintf(stderr, "#Synthetic Dataset\n");
        fprintf(stderr, "Number of Locations: %d\n", N);
        fprintf(stderr, "#Threads per node: %d\n", ncores);
        fprintf(stderr, "#GPUs: %d\n", gpus);
        if (precision == 1)
            fprintf(stderr, "#Double Precision!\n");
        else if (precision == 0)
            fprintf(stderr, "#Single Precision!\n");
        else if (precision == 2)
            fprintf(stderr, "#Single/Double Precision!\n");
        fprintf(stderr, "#Dense Tile Size: %d\n", ts);
        fprintf(stderr, "#exact computation\n");
        fprintf(stderr, "p=%d, q=%d\n", p_grid, q_grid);
        fprintf(stderr, "***************************************************\n");
#if defined(CHAMELEON_USE_MPI)
        }
#endif
        is_printed = true;
    }
}

#endif //EXAGEOSTATCPP_UTILS_HPP
