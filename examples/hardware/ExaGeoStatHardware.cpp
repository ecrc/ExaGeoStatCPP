
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file TestExaGeoStatHardware.cpp
 * @briefExample for the ExaGeoStatHardware class in the ExaGeoStat software package.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @date 2024-02-14
**/

#include <hardware/ExaGeoStatHardware.hpp>
#include <utilities/Logger.hpp>

using namespace exageostat::configurations;

/**
 * @brief Main entry point for demonstrating the ExaGeoStatHardware class usage.
 * @details This program demonstrates how to initialize the ExaGeoStatHardware class using configuration settings derived from command-line arguments.
 * It showcases the initialization process for hardware configurations, including the computation mode, number of cores, and number of GPUs,
 * highlighting the simplicity and effectiveness of managing hardware resources in ExaGeoStat.
 * @param argc Number of command-line arguments.
 * @param argv Array of command-line argument strings.
 * @return An integer indicating the success or failure of the program. A return value of 0 indicates success, while any non-zero value indicates failure.
 *
 */
int main(int argc, char **argv) {

    LOGGER("** Example of Hardware **")

    // Initialize Configuration
    Configurations configuration;
    configuration.InitializeArguments(argc, argv);

    // Initialize Hardware
    auto hardware = ExaGeoStatHardware(configuration.GetComputation(), configuration.GetCoresNumber(),
                                       configuration.GetGPUsNumbers());

    return 0;
}