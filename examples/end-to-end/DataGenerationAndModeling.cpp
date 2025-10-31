
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file DataGenerationAndModeling.cpp
 * @brief This program either generates synthetic data using the ExaGeoStat library, or reads an CSV file.
 * @details The program takes command line arguments to configure the data generation.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @date 2024-02-04
**/

#include <api/ExaGeoStat.hpp>

using namespace exageostat::api;
using namespace exageostat::configurations;

/**
 * @brief Main entry point for the Data Generation & Data Modeling program.
 * @details This function either generates synthetic data using the ExaGeoStat library, or reads an CSV file, then models the loaded data.
 * @param[in] argc The number of command line arguments.
 * @param[in] argv An array of command line argument strings.
 * @return An integer indicating the success or failure of the program. A return value of 0 indicates success, while any non-zero value indicates failure.
 *
 */
int main(int argc, char **argv) {

    // Create a new configurations object.
    Configurations configurations;
    //  Initialize the arguments with the provided command line arguments
    configurations.InitializeArguments(argc, argv);
    // Initialize the ExaGeoStat Hardware based on runtime type
#if DEFAULT_RUNTIME
    // StarPU/CHAMELEON mode
    auto hardware = ExaGeoStatHardware(configurations.GetComputation(), configurations.GetCoresNumber(),
                                       configurations.GetGPUsNumbers(), configurations.GetPGrid(),
                                       configurations.GetQGrid());
#else
    // PaRSEC mode
    auto hardware = ExaGeoStatHardware(configurations);
#endif
    // Load data by either read from file or create synthetic data.
    std::unique_ptr<ExaGeoStatData<double>> data;
    ExaGeoStat<double>::ExaGeoStatLoadData(configurations, data);
    // Modeling module.
    ExaGeoStat<double>::ExaGeoStatDataModeling(configurations, data);

    return 0;
}

