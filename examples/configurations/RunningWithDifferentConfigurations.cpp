
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file RunningWithDifferentConfigurations.cpp
 * @brief Demonstrates running ExaGeoStat with various configurations.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @date 2024-01-04
**/

#include <api/ExaGeoStat.hpp>
#include <utilities/Logger.hpp>

using namespace exageostat::api;
using namespace exageostat::configurations;

/**
 * @brief Main entry point for the Data Generation & Data Modeling program.
 * @details Initializes configuration settings, loads data, and performs data modeling with ExaGeoStat.
 * The program demonstrates running ExaGeoStat with two different sets of configurations to showcase the software's capability in handling various statistical models and computational settings.
 * @return An integer indicating the success or failure of the program. A return value of 0 indicates success, while any non-zero value indicates failure.
 *
 */
int main() {

    // Create a new configurations object.
    Configurations configurations;
    int N = 16;
    configurations.SetProblemSize(N);
    configurations.SetKernelName("UnivariateMaternStationary");
    configurations.SetDenseTileSize(9);
    configurations.SetMaxMleIterations(10);
    configurations.SetTolerance(4);
    std::vector<double> lb{0.1, 0.1, 0.1};
    configurations.SetLowerBounds(lb);

    std::vector<double> ub{5, 5, 5};
    configurations.SetUpperBounds(ub);

    std::vector<double> initial_theta{1, 0.1, 0.5};
    configurations.SetInitialTheta(initial_theta);

    // Initialize the ExaGeoStat Hardware
    ExaGeoStatHardware hardware(configurations.GetComputation(), configurations.GetCoresNumber(),configurations.GetGPUsNumbers());
    // Load data by either read from file or create synthetic data.
    std::unique_ptr<ExaGeoStatData<double>> data;
    ExaGeoStat<double>::ExaGeoStatLoadData(configurations, data);
    // Modeling module.
    ExaGeoStat<double>::ExaGeoStatDataModeling(configurations, data);

    LOGGER("")
    LOGGER("ANOTHER CONFIGURATIONS\n")
    // Create a new configurations object.
    Configurations configurations2;
    N = 10;
    configurations2.SetProblemSize(N);
    configurations2.SetKernelName("UnivariateMaternStationary");
    configurations2.SetDenseTileSize(5);
    configurations2.SetMaxMleIterations(2);
    configurations2.SetTolerance(3);
    configurations2.SetLowerBounds(lb);
    configurations2.SetUpperBounds(ub);
    configurations2.SetInitialTheta(initial_theta);

    // Load data by either read from file or create synthetic data.
    ExaGeoStat<double>::ExaGeoStatLoadData(configurations2, data);
    // Modeling module.
    ExaGeoStat<double>::ExaGeoStatDataModeling(configurations2, data);

    return 0;
}
