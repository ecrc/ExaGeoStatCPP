
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file SyntheticLocationsGeneration.cpp
 * @brief This file contains the main function for generating synthetic Locations for ExaGeoStat
 * @version 1.0.0
 * @author Mahmoud ElKarargy
 * @date 2023-03-04
**/

#include <iostream>

#include <data-generators/DataGenerator.hpp>
#include <api/ExaGeoStat.hpp>

using namespace std;

using namespace exageostat::configurations;
using namespace exageostat::generators;
using namespace exageostat::common;
using namespace exageostat::hardware;

/**
 * @brief The main function of the program.
 * @details This function generates synthetic data for ExaGeoStat using the provided command line arguments.
 * @param[in] argc The number of command line arguments.
 * @param[in] argv The command line arguments.
 * @return The status code of the program.
 */

int main(int argc, char **argv) {

    // Create a new synthetic_data_configurations object with the provided command line arguments
    Configurations synthetic_data_configurations;
    synthetic_data_configurations.InitializeArguments(argc, argv);

    // initialize ExaGeoStat Hardware.
    auto hardware = ExaGeoStatHardware(synthetic_data_configurations.GetComputation(),
                                       synthetic_data_configurations.GetCoresNumber(),
                                       synthetic_data_configurations.GetGPUsNumbers());

    // Create a unique pointer to a DataGenerator object
    unique_ptr<DataGenerator<double>> synthetic_generator = DataGenerator<double>::CreateGenerator(
            synthetic_data_configurations);

    // Initialize the locations of the generated data
    auto *locations = synthetic_generator->CreateLocationsData(synthetic_data_configurations);

    // Define a struct to hold pointers to the x, y, and z coordinates of the generated data
    struct DataPointers {
        double *x;
        double *y;
        double *z;
    } data_pointers{};

    // Set the pointers in the DataPointers struct to the location coordinates of the generated data
    data_pointers.x = locations->GetLocationX();
    data_pointers.y = locations->GetLocationY();
    data_pointers.z = locations->GetLocationZ();

    // Print the generated location coordinates
    cout << "Generated Locations are .. " << endl;
    int timeSlot;
    if (synthetic_data_configurations.GetDimension() != DimensionST) {
        timeSlot = 1;
    } else {
        timeSlot = synthetic_data_configurations.GetTimeSlot();
    }
    for (auto i = 0; i < synthetic_data_configurations.GetProblemSize() * timeSlot; i++) {
        printf("X: %.18f Y: %.18f", data_pointers.x[i], data_pointers.y[i]);
        if (synthetic_data_configurations.GetDimension() != Dimension2D) {
            cout << " Z: " << data_pointers.z[i];
        }
        printf("\n");
    }

    return 0;
}