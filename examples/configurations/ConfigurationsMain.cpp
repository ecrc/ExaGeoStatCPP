
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
* @file ConfigurationsMain.cpp
* @brief Demonstrates how to use the Configurations class from the ExaGeoStat software package.
* @details This file demonstrates how to use the Configurations class from the ExaGeoStat software package
* to obtain user-defined configurations for generating synthetic data.
* @version 1.0.0
* @author Sameh Abdulah
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

    // Obtain user-defined configurations and print them to the console.
    int N = configurations.GetProblemSize();
    if (N != 0) {
        cout << "You set N by: " << N << endl;
    }

    string kernel = configurations.GetKernelName();
    if (!kernel.empty()) {
        cout << "You set Kernel by: " << kernel << endl;
    }

    Dimension dimension = configurations.GetDimension();
    if (dimension == Dimension2D) {
        cout << "You set Dimension by: 2D" << endl;
    } else if (dimension == Dimension3D) {
        cout << "You set Dimension by: 3D" << endl;
    } else if (dimension == DimensionST) {
        cout << "You set Dimension by: ST" << endl;
    }

    int p_grid = configurations.GetPGrid();
    if (p_grid != 0) {
        cout << "You set P by: " << p_grid << endl;
    }

    int time_slot = configurations.GetTimeSlot();
    if (time_slot != 0) {
        cout << "You set time slot by: " << time_slot << endl;
    }

    Computation computation = configurations.GetComputation();
    if (computation == EXACT_DENSE) {
        cout << "You set Computation to: EXACT" << endl;
    } else if (computation == DIAGONAL_APPROX) {
        cout << "You set Computation to: DIAGONAL APPROX" << endl;
    } else if (computation == TILE_LOW_RANK) {
        cout << "You set Computation to: TILE LOW RANK" << endl;
    }

    Precision precision = configurations.GetPrecision();
    if (precision == SINGLE) {
        cout << "You set precision to: SINGLE" << endl;
    } else if (precision == DOUBLE) {
        cout << "You set precision to: DOUBLE" << endl;
    } else if (precision == MIXED) {
        cout << "You set precision to: MIXED PRECISION" << endl;
    }

    int seed = configurations.GetSeed();
    cout << "You set seed to: " << seed << endl;

    int run_mode = Configurations::GetRunMode();
    if (run_mode == RunMode::VERBOSE_MODE) {
        cout << " You set run mode to VERBOSE." << endl;
    } else {
        cout << "You set run mode to STANDARD." << endl;
    }

    VERBOSE("VERBOSE ACTIVATED")

    return 0;
}