// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file Utils.hpp
 * @version 1.0.0
 * @brief This file contains common functions used in ExaGeoStat software package.
 * @details These functions include so far the VERBOSE macro that prints a message
 * to the console if the verbosity setting is set to "verbose mode.
 * @author Sameh Abdulah
 * @date 2023-03-21
**/

#ifndef EXAGEOSTATCPP_UTILS_HPP
#define EXAGEOSTATCPP_UTILS_HPP

#include <iostream>
#include <string>

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
        std::cout << msg << std::endl; \
        exit(EXIT_FAILURE);\
    }

#endif //EXAGEOSTATCPP_UTILS_HPP
