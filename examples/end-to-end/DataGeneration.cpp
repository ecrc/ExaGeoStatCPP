
/*
 * Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
 * All rights reserved.
 * ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).
 */

/**
 * @file DataGeneration.cpp
 * @brief This program generates synthetic data using the ExaGeoStat library.
 * @details The program takes command line arguments to configure the data generation.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-05-30
**/

#include <configurations/data-generation/concrete/SyntheticDataConfigurations.hpp>
#include <api/ExaGeoStat.hpp>

using namespace exageostat::configurations::data_configurations;
using namespace exageostat::api;

/**
 * @brief Main entry point for the DataGeneration program.
 * @details This function generates synthetic data using the ExaGeoStat library.
 * @param argc The number of command line arguments.
 * @param argv An array of command line argument strings.
 * @return An integer indicating the success or failure of the program.
 */
int main(int argc, char **argv) {

    // Create a new synthetic_data_configurations object with the provided command line arguments
    SyntheticDataConfigurations synthetic_data_configurations(argc, argv);

    std::cout << "** Initialise ExaGeoStat hardware ** " << std::endl;
    // Initialise ExaGeoStat hardware with the selected number of cores and  gpus.
    ExaGeoStat<double>::ExaGeoStatInitializeHardware(&synthetic_data_configurations);

    std::cout << "** Generate ExaGeoStat data ** " << std::endl;
    ExaGeoStat<double>::ExaGeoStatGenerateData(&synthetic_data_configurations);

    std::cout << "** Finalize ExaGeoStat hardware ** " << std::endl;
    // Finalise ExaGeoStat context.
    ExaGeoStat<double>::ExaGeoStatFinalizeHardware(&synthetic_data_configurations);

    return 0;
}