
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file Utils.hpp
 * @brief This file contains common functions used in ExaGeoStat software package.
 * @details These functions include so far the VERBOSE macro that prints a message
 * to the console if the verbosity setting is set to "verbose mode.
 * @version 1.0.0
 * @author Mahmoud ElKarargy
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
 * Timing macro to start timing.
 */
#define START_TIMING(t) auto t##_start = std::chrono::high_resolution_clock::now()

/**
 * Timing macro to stop timing.
 */
#define STOP_TIMING(t) auto t##_end = std::chrono::high_resolution_clock::now(); \
                    t = std::chrono::duration_cast<std::chrono::duration<double>>(t##_end - t##_start).count()

#endif //EXAGEOSTATCPP_UTILS_HPP
