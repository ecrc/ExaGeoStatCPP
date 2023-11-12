
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
#include <sys/time.h>

#include <common/Definitions.hpp>
#include <configurations/Configurations.hpp>
#include <helpers/CommunicatorMPI.hpp>

/**
 * DEFAULT_PRECISION the value of the default C++ std::cout number of precision.
 */
#define DEFAULT_PRECISION 6

/**
* Verbose macro for logging and debugging mode
*/
#define VERBOSE(msg) \
    if(exageostat::configurations::Configurations::GetVerbosity() == exageostat::common::Verbose::DETAILED_MODE && exageostat::helpers::CommunicatorMPI::GetInstance()->GetRank()) \
        std::cout << "\t\t\t " << msg << std::endl;

/**
 * LOGGER_1 macro for logging outputs with double taps and new line at the end.
 */

#define LOGGER_1(msg) \
    if(!(exageostat::configurations::Configurations::GetVerbosity() == exageostat::common::Verbose::QUIET_MODE) && exageostat::helpers::CommunicatorMPI::GetInstance()->GetRank()) \
        std::cout << "\t\t " << std::fixed << std::setprecision(DEFAULT_PRECISION) << msg << std::endl;

/**
 * LOGGER_2 macro for logging outputs with double taps and without new line at the end.
 */
#define LOGGER_2(msg, A) \
    if(!(exageostat::configurations::Configurations::GetVerbosity() == exageostat::common::Verbose::QUIET_MODE) && exageostat::helpers::CommunicatorMPI::GetInstance()->GetRank()) \
        std::cout << "\t\t " << std::fixed << std::setprecision(DEFAULT_PRECISION) << msg;

/**
 * LOGGER_CONTROL is The internal macro that simply strips the excess and ends up with the required macro
 */
#define LOGGER_CONTROL(x, A, B, FUNC, ...)  FUNC

/**
 * LOGGER macro that's called, Used to logging outputs
 */
#define LOGGER(...)    LOGGER_CONTROL(,##__VA_ARGS__, \
                       LOGGER_2(__VA_ARGS__),         \
                       LOGGER_1(__VA_ARGS__),         \
                                         )

/**
* LOGGER_PRECISION_1 macro for logging outputs without any taps, without new line at the end and with customized precision.
*/
#define LOGGER_PRECISION_1(msg, precision) \
    if(!(exageostat::configurations::Configurations::GetVerbosity() == exageostat::common::Verbose::QUIET_MODE) && exageostat::helpers::CommunicatorMPI::GetInstance()->GetRank()) \
        std::cout << std::fixed << std::setprecision(precision) << msg;

/**
* LOGGER_PRECISION macro for logging outputs without any taps, without new line at the end and with default C++ precision.
*/
#define LOGGER_PRECISION_2(msg) \
    if(!(exageostat::configurations::Configurations::GetVerbosity() == exageostat::common::Verbose::QUIET_MODE) && exageostat::helpers::CommunicatorMPI::GetInstance()->GetRank()) \
        std::cout << std::fixed << std::setprecision(DEFAULT_PRECISION) << msg;

/**
 * LOGGER_CONTROL is The internal macro that simply strips the excess and ends up with the required macro
 */
#define LOGGER_PRECISION_CONTROL(x, A, B, FUNC, ...)  FUNC

/**
 * LOGGER macro that's called, Used to logging outputs
 */
#define LOGGER_PRECISION(...)    LOGGER_PRECISION_CONTROL(,##__VA_ARGS__, \
                                 LOGGER_PRECISION_1(__VA_ARGS__),         \
                                 LOGGER_PRECISION_2(__VA_ARGS__),         \
                                                 )

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
