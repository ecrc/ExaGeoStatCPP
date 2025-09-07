
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file StageZeroDataGenerator.cpp
 * @brief Example to run Stage Zero data generation (mean-trend pipeline).
 * @version 2.0.0
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
    
    // Initialize the ExaGeoStat Hardware based on runtime type
#if DEFAULT_RUNTIME
    // StarPU/CHAMELEON mode - use constructor that initializes CHAMELEON
    auto hardware = ExaGeoStatHardware(configurations.GetComputation(), 
                                       configurations.GetCoresNumber(),
                                       configurations.GetGPUsNumbers(), 
                                       configurations.GetPGrid(), 
                                       configurations.GetQGrid());
#else
    // PaRSEC mode - use constructor that initializes PaRSEC
    auto hardware = ExaGeoStatHardware(configurations);
#endif
    
    std::unique_ptr<ExaGeoStatData<double>> data;
    // Generate Stage Zero mean-trend data
    ExaGeoStat<double>::ExaGeoStatGenerateMeanTrendData(configurations, data);

    // Finalize Hardware
    hardware.FinalizeHardware();

    return 0;
}
