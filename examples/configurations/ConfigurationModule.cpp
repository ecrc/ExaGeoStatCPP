
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
* @file ConfigurationsMain.cpp
* @brief Demonstrates how to use the Configurations class from the ExaGeoStat software package.
* @details This file demonstrates how to use the Configurations class from the ExaGeoStat software package
* to obtain user-defined configurations for generating synthetic data.
* @version 1.0.0
* @author Mahmoud ElKarargy
* @date 2023-01-31
*
**/

#include <iostream>

#include <common/Utils.hpp>
#include <configurations/Configurations.hpp>

using namespace std;

using namespace exageostat::common;
using namespace exageostat::configurations;

/**
 * @brief The main function of the program.
 *
 * This function demonstrates how to use the SyntheticDataConfigurations class from the ExaGeoStat software package
 * to obtain user-defined configurations for generating synthetic data.
 *
 * @param[in] argc The number of command line arguments.
 * @param[in] argv The command line arguments.
 * @return The status code of the program.
 */
int main(int argc, char **argv) {

    // Create an instance of the SyntheticDataConfigurations class with user-defined configurations.
    Configurations configurations;
    configurations.InitializeArguments(argc, argv);

    LOGGER("** These are some examples of the common arguments needed between all modules of ExaGeoStat **")

    // Obtain user-defined configurations and print them to the console.
    int n = configurations.GetProblemSize();
    if (n != 0) {
        LOGGER("You set N by: " << n)
    }

    string kernel = configurations.GetKernelName();
    if (!kernel.empty()) {
        LOGGER("You set Kernel by: " << kernel)
    }

    Dimension dimension = configurations.GetDimension();
    if (dimension == Dimension2D) {
        LOGGER("You set Dimension by: 2D")
    } else if (dimension == Dimension3D) {
        LOGGER("You set Dimension by: 3D")
    } else if (dimension == DimensionST) {
        LOGGER("You set Dimension by: ST")
    }

    int p_grid = configurations.GetPGrid();
    if (p_grid != 0) {
        LOGGER("You set P by: " << p_grid)
    }

    int time_slot = configurations.GetTimeSlot();
    if (time_slot != 0) {
        LOGGER("You set time slot by: " << time_slot)
    }

    Computation computation = configurations.GetComputation();
    if (computation == EXACT_DENSE) {
        LOGGER("You set Computation to: EXACT")
    } else if (computation == DIAGONAL_APPROX) {
        LOGGER("You set Computation to: DIAGONAL APPROX")
    } else if (computation == TILE_LOW_RANK) {
        LOGGER("You set Computation to: TILE LOW RANK")
    }

    Precision precision = configurations.GetPrecision();
    if (precision == SINGLE) {
        LOGGER("You set precision to: SINGLE")
    } else if (precision == DOUBLE) {
        LOGGER("You set precision to: DOUBLE")
    } else if (precision == MIXED) {
        LOGGER("You set precision to: MIXED PRECISION")
    }

    int seed = configurations.GetSeed();
    LOGGER("You set seed to: " << seed)

    int run_mode = Configurations::GetVerbosity();
    if (run_mode == Verbose::DETAILED_MODE) {
        LOGGER(" You set run mode to VERBOSE.")
    } else {
        LOGGER("You set run mode to STANDARD.")
    }

    VERBOSE("VERBOSE ACTIVATED")
    return 0;
}