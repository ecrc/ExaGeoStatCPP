
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file DataGeneration.cpp
 * @brief This program either generates synthetic data using the ExaGeoStat library, or reads an CSV file containing real data.
 * @details The program takes command line arguments to configure the data generation.
 * @version 1.0.0
 * @author Mahmoud ElKarargy
 * @date 2023-05-30
**/

#include <common/Utils.hpp>
#include <api/ExaGeoStat.hpp>

using namespace exageostat::configurations;
using namespace exageostat::api;
using namespace exageostat::hardware;
using namespace exageostat::dataunits;

/**
 * @brief Main entry point for the DataGeneration program.
 * @details This function either generates synthetic data using the ExaGeoStat library, or reads an CSV file containing real data.
 * @param[in] argc The number of command line arguments.
 * @param[in] argv An array of command line argument strings.
 * @return An integer indicating the success or failure of the program.
 */
int main(int argc, char **argv) {

    // Create a new configurations object.
    Configurations configurations;
    //  Initialize the arguments with the provided command line arguments
    configurations.InitializeArguments(argc, argv);
    LOGGER("** initialize ExaGeoStat hardware ** ")
    auto hardware = ExaGeoStatHardware(configurations.GetComputation(), configurations.GetCoresNumber(),
                                       configurations.GetGPUsNumbers());
    LOGGER("** Create ExaGeoStat data **")
    ExaGeoStatData<double> data;
    LOGGER("** Generate ExaGeoStat data ** ")
    ExaGeoStat<double>::ExaGeoStatLoadData(hardware, configurations, data);
    LOGGER("** All example stages have been completed successfully ** ")

    return 0;
}