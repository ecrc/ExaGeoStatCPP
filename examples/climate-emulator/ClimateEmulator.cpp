
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ClimateEmulator.cpp
 * @brief example of climate emulator.
 * @version 2.0.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2024-09-23
**/

#include <configurations/Configurations.hpp>
#include <hardware/ExaGeoStatHardware.hpp>
#include <api/ExaGeoStat.hpp>

using namespace exageostat::configurations;
using namespace exageostat::api;

int main(int argc, char **argv) {

    // Create a new configurations object.
    Configurations configurations;
    // Initialize the arguments with the provided command line arguments
    configurations.InitializeArguments(argc, argv);
    // Initialize the ExaGeoStat Hardware
    auto hardware = ExaGeoStatHardware(configurations);
    // Create a unique pointer to hold the data.
    std::unique_ptr<ExaGeoStatData<double>> data;
    // Load the data, either by reading from a file or generating synthetic data.
    ExaGeoStat<double>::ExaGeoStatLoadData(configurations, data);
    // Perform data modeling.
    ExaGeoStat<double>::ExaGeoStatDataModeling(configurations, data);

    return 0;
}
