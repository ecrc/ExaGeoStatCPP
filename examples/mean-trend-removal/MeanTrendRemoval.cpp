
// Copyright (c) 2017-2025 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file MeanTrendRemovalDataGenerator.cpp
 * @brief Example to run Mean Trend Removal data generation (mean-trend pipeline).
 * @version 2.0.0
 * @author Ali Hakam
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2025-11-12
 **/

#include <configurations/Configurations.hpp>
#include <hardware/ExaGeoStatHardware.hpp>
#include <api/ExaGeoStat.hpp>

using namespace exageostat::configurations;
using namespace exageostat::api;

int main(int argc, char **argv) {

    // Create a configurations object.
    Configurations configurations;
    // Initialize the arguments with the provided command line arguments
    configurations.InitializeArguments(argc, argv);
    
    // Initialize the ExaGeoStat Hardware
    auto hardware = ExaGeoStatHardware(configurations);
    
    std::unique_ptr<ExaGeoStatData<double>> data;
    // Generate Mean Trend Removal mean-trend data
    ExaGeoStat<double>::ExaGeoStatGenerateMeanTrendData(configurations, data);

    // Finalize Hardware
    hardware.FinalizeHardware();

    return 0;
}
