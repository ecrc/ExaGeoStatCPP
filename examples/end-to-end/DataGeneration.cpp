
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

#include <api/ExaGeoStat.hpp>
#include <configurations/Configurations.hpp>

using namespace exageostat::configurations;
using namespace exageostat::api;

/**
 * @brief Main entry point for the DataGeneration program.
 * @details This function generates synthetic data using the ExaGeoStat library.
 * @param argc The number of command line arguments.
 * @param argv An array of command line argument strings.
 * @return An integer indicating the success or failure of the program.
 */
int main(int argc, char **argv) {

//    // Create a new configurations object with the provided command line arguments
//    auto conf = Configurations::GetConfigurations()->InitializeArguments(argc, argv);
//
//    std::cout << "** Initialise ExaGeoStat hardware ** " << std::endl;
//    // Initialise ExaGeoStat hardware with the selected number of cores and  gpus.
//    ExaGeoStat<double>::ExaGeoStatInitializeHardware();
//
//        // class ExageostatData (descriptors that will be generated)
//        //class LocationData
//        //struct {Exa Loc } generation_result
//    std::cout << "** Generate ExaGeoStat data ** " << std::endl;
//    auto exageostat_data = ExaGeoStat<double>::ExaGeoStatGenerateData(conf);
//
//    std::cout << "** Finalize ExaGeoStat hardware ** " << std::endl;
//    ExaGeoStat<double>::ExaGeoStatFinalizeHardware();

    return 0;
}