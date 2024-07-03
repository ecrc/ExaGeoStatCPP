
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file Logger.hpp
 * @brief Provides logging and timing macros for debugging and profiling.
 * @details Defines macros for verbose logging, various levels of logging, and timing.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2024-02-04
**/

#ifndef EXAGEOSTATCPP_LOGGER_HPP
#define EXAGEOSTATCPP_LOGGER_HPP

#include <iostream>
#include <string>
#include <sys/time.h>

#include <common/Definitions.hpp>
#include <configurations/Configurations.hpp>
#include <helpers/CommunicatorMPI.hpp>

/**
 * @def DEFAULT_PRECISION
 * @brief The value of the default C++ std::std::cout number of precision.
 */
#define DEFAULT_PRECISION 6

/**
 * @def VERBOSE_PRECISION
 * @brief A fixed precision for the verbose mode.
 */
#define VERBOSE_PRECISION 12

/**
 * @def VERBOSE_1(msg)
 * @brief VERBOSE_1 macro for logging outputs in verbose mode, with double taps at the beginning, a fixed precision and a new line at the end.
 */
#define VERBOSE_1(msg) \
    if(exageostat::configurations::Configurations::GetVerbosity() == exageostat::common::Verbose::DETAILED_MODE && \
      !exageostat::helpers::CommunicatorMPI::GetInstance()->GetRank()) { \
        std::ostringstream oss; \
        oss << "\t\t " << std::fixed << std::setprecision(VERBOSE_PRECISION) << msg << std::endl; \
        std::cout << oss.str(); \
    }

/**
 * @def VERBOSE_2(msg, A)
 * @brief VERBOSE_2 macro for logging outputs in verbose mode, with double taps at the beginning, a fixed precision and without a new line at the end.
 */
#define VERBOSE_2(msg, A) \
    if(exageostat::configurations::Configurations::GetVerbosity() == exageostat::common::Verbose::DETAILED_MODE && \
       !exageostat::helpers::CommunicatorMPI::GetInstance()->GetRank()) { \
        std::ostringstream oss; \
        oss << "\t\t " << std::fixed << std::setprecision(VERBOSE_PRECISION) << msg ; \
        std::cout << oss.str(); \
    }

/**
 * @def VERBOSE_CONTROL(x, A, B, FUNC, ...)
 * @brief VERBOSE_CONTROL is The internal macro that simply strips the excess and ends up with the required macro
 */
#define VERBOSE_CONTROL(x, A, B, FUNC, ...) FUNC

/**
 * @def VERBOSE(...)
 * @brief VERBOSE macro that's called, Used to logging outputs in a verbose mode.
 */
#define VERBOSE(...) VERBOSE_CONTROL(,##__VA_ARGS__, \
                       VERBOSE_2(__VA_ARGS__),         \
                       VERBOSE_1(__VA_ARGS__),       \
                                         )
/**
 * @def LOGGER_1(msg)
 * @brief LOGGER_1 macro for logging outputs with double taps and new line at the end.
 */
#define LOGGER_1(msg) \
    if(!(exageostat::configurations::Configurations::GetVerbosity() == exageostat::common::Verbose::QUIET_MODE) && !exageostat::helpers::CommunicatorMPI::GetInstance()->GetRank()){ \
        std::ostringstream oss; \
        oss << "\t\t " << std::fixed << std::setprecision(DEFAULT_PRECISION) << msg << std::endl; \
        std::cout << oss.str();                                                                                                                                        \
    }

/**
 * @def LOGGER_2(msg, A)
 * @brief LOGGER_2 macro for logging outputs with double taps and without new line at the end.
 */
#define LOGGER_2(msg, A) \
    if(!(exageostat::configurations::Configurations::GetVerbosity() == exageostat::common::Verbose::QUIET_MODE) && !exageostat::helpers::CommunicatorMPI::GetInstance()->GetRank()){ \
        std::ostringstream oss;                                                                                                                                                     \
        oss << "\t\t " << std::fixed << std::setprecision(DEFAULT_PRECISION) << msg;                 \
        std::cout << oss.str(); \
}
/**
 * @def LOGGER_CONTROL(x, A, B, FUNC, ...)
 * @brief LOGGER_CONTROL is The internal macro that simply strips the excess and ends up with the required macro
 */
#define LOGGER_CONTROL(x, A, B, FUNC, ...)  FUNC

/**
 * @def LOGGER(...)
 * @brief LOGGER macro that's called, Used to logging outputs.
 */
#define LOGGER(...)    LOGGER_CONTROL(,##__VA_ARGS__, \
                       LOGGER_2(__VA_ARGS__),         \
                       LOGGER_1(__VA_ARGS__),         \
                                         )

/**
 * @def LOGGER_PRECISION_1(msg, precision)
 * @brief LOGGER_PRECISION_1 macro for logging outputs without any taps, without new line at the end, and with customized precision.
 */
#define LOGGER_PRECISION_1(msg, precision) \
    if(!(exageostat::configurations::Configurations::GetVerbosity() == exageostat::common::Verbose::QUIET_MODE) && !exageostat::helpers::CommunicatorMPI::GetInstance()->GetRank()){ \
        std::ostringstream oss;                                                                                                                                                     \
        oss << std::fixed << std::setprecision(precision) << msg;                                   \
        std::cout << oss.str(); \
    }

/**
 * @def LOGGER_PRECISION_2(msg)
 * @brief LOGGER_PRECISION_2 macro for logging outputs without any taps, without new line at the end, and with default C++ precision.
 */
#define LOGGER_PRECISION_2(msg) \
    if(!(exageostat::configurations::Configurations::GetVerbosity() == exageostat::common::Verbose::QUIET_MODE) && !exageostat::helpers::CommunicatorMPI::GetInstance()->GetRank()) {\
        std::ostringstream oss;                                                                                                                                                     \
        oss << std::fixed << std::setprecision(DEFAULT_PRECISION) << msg;                                   \
        std::cout << oss.str(); \
    }
/**
 * @def LOGGER_PRECISION_CONTROL
 * @brief is The internal macro that simply strips the excess and ends up with the required macro
 */
#define LOGGER_PRECISION_CONTROL(x, A, B, FUNC, ...)  FUNC

/**
 * @def LOGGER_PRECISION(...)
 * @brief LOGGER_PRECISION macro that's called, Used for logging outputs with precision.
 */
#define LOGGER_PRECISION(...)    LOGGER_PRECISION_CONTROL(,##__VA_ARGS__, \
                                 LOGGER_PRECISION_1(__VA_ARGS__),         \
                                 LOGGER_PRECISION_2(__VA_ARGS__),         \
                                                 )

/**
 * @def START_TIMING(t)
 * @brief Timing macro to start timing.
 */
#define START_TIMING(t) auto t##_start = std::chrono::high_resolution_clock::now()

/**
 * @def STOP_TIMING(t)
 * @brief Timing macro to stop timing.
 */
#define STOP_TIMING(t) auto t##_end = std::chrono::high_resolution_clock::now(); \
                    t = std::chrono::duration_cast<std::chrono::duration<double>>(t##_end - t##_start).count()

#endif //EXAGEOSTATCPP_LOGGER_HPP
