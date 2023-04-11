/*
 * Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
 * Copyright (C) 2023 by Brightskies inc,
 * All rights reserved.
 * ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).
 */

/**
* @file main_synthetic.cpp
* @brief Demonstrates how to use the SyntheticDataConfigurations class from the ExaGeoStat software package.
* @version 1.0.0
* @author Sameh Abdulah
* @date 2023-01-31
*
* This file demonstrates how to use the SyntheticDataConfigurations class from the ExaGeoStat software package
* to obtain user-defined configurations for generating synthetic data.
*
**/

#include <iostream>
#include <configurations/data-generation/concrete/SyntheticDataConfigurations.hpp>

using namespace std;
using namespace exageostat::configurations::data_configurations;
using namespace exageostat::common;

/**
 * @brief The main function of the program.
 *
 * This function demonstrates how to use the SyntheticDataConfigurations class from the ExaGeoStat software package
 * to obtain user-defined configurations for generating synthetic data.
 *
 * @param argc The number of command line arguments.
 * @param argv The command line arguments.
 * @return The status code of the program.
 */
int main(int argc, char **argv) {

    // Create an instance of the SyntheticDataConfigurations class with user-defined configurations.
    auto synthetic_data_configurations = new SyntheticDataConfigurations(argc, argv);

    // Obtain user-defined configurations and print them to the console.
    int N = synthetic_data_configurations->GetProblemSize();
    if (N != 0) {
        cout << "You set N by: " << N << endl;
    }

    string kernel = synthetic_data_configurations->GetKernel();
    if (!kernel.empty()) {
        cout << "You set Kernel by: " << kernel << endl;
    }

    Dimension dimension = synthetic_data_configurations->GetDimension();
    if (dimension == Dimension2D) {
        cout << "You set Dimension by: 2D" << endl;
    } else if (dimension == Dimension3D) {
        cout << "You set Dimension by: 3D" << endl;
    } else if (dimension == DimensionST) {
        cout << "You set Dimension by: ST" << endl;
    }

    int p_grid = synthetic_data_configurations->GetPGrid();
    if (p_grid != 0) {
        cout << "You set P by: " << p_grid << endl;
    }

    int time_slot = synthetic_data_configurations->GetTimeSlot();
    if (time_slot != 0) {
        cout << "You set time slot by: " << time_slot << endl;
    }

    Computation computation = synthetic_data_configurations->GetComputation();
    if (computation == EXACT_DENSE) {
        cout << "You set Computation to: EXACT" << endl;
    } else if (computation == DIAGONAL_APPROX) {
        cout << "You set Computation to: DIAGONAL APPROX" << endl;
    } else if (computation == TILE_LOW_RANK) {
        cout << "You set Computation to: TILE LOW RANK" << endl;
    }

    Precision precision = synthetic_data_configurations->GetPrecision();
    if (precision == SINGLE) {
        cout << "You set precision to: SINGLE" << endl;
    } else if (precision == DOUBLE) {
        cout << "You set precision to: DOUBLE" << endl;
    } else if (precision == MIXED) {
        cout << "You set precision to: MIXED PRECISION" << endl;
    }
    // Clean up memory
    delete synthetic_data_configurations;
    return 0;
}