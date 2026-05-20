
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file DataModeling.cpp
 * @brief This program models data using the ExaGeoStat library.
 * @details The program takes command line arguments and example variables to configure the data modeling.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @date 2024-02-04
**/

#include <api/ExaGeoStat.hpp>

using namespace exageostat::api;
using namespace exageostat::configurations;

/**
 * @brief Main entry point for the Data Modeling program.
 * @details This example illustrates the process of data modeling using the ExaGeoStat library's CHAMELEON descriptor framework.
 * It involves configuring parameters such as problem size and computation mode, initializing hardware resources, setting up matrices for descriptors,
 * and creating location information. The ExaGeoStatDataModeling function is then called to perform geo statistical analysis. The example showcases
 * the library's efficiency in handling large spatial datasets while efficiently utilizing hardware resources..
 * @param[in] argc The number of command line arguments.
 * @param[in] argv An array of command line argument strings.
 * @return An integer indicating the success or failure of the program. A return value of 0 indicates success, while any non-zero value indicates failure.
 *
 */
int main(int argc, char **argv) {
    // Create a new data_modeling_configurations object with the provided command line arguments and example variables
    Configurations configurations;
    configurations.InitializeArguments(argc, argv);
    // initialize ExaGeoStat hardware with the selected number of cores and  gpus.
#if DEFAULT_RUNTIME
    // StarPU/CHAMELEON mode
    auto hardware = ExaGeoStatHardware(configurations.GetComputation(), configurations.GetCoresNumber(),
                                       configurations.GetGPUsNumbers(), configurations.GetPGrid(),
                                       configurations.GetQGrid());
#else
    // PaRSEC mode
    auto hardware = ExaGeoStatHardware(configurations);
#endif

    // Load data from CSV file (via --data_path) or generate synthetic data.
    std::unique_ptr<ExaGeoStatData<double>> data;
    ExaGeoStat<double>::ExaGeoStatLoadData(configurations, data);

    // Modeling module.
    ExaGeoStat<double>::ExaGeoStatDataModeling(configurations, data);

    return 0;
}
