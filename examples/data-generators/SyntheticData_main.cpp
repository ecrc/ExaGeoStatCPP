
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file synthetic_generator.cpp
 * @brief This file contains the main function for generating synthetic data for ExaGeoStat.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-03-04
**/

#include <iostream>
#include <data-generators/DataGenerator.hpp>

using namespace exageostat::generators;
using namespace exageostat::dataunits;
using namespace exageostat::common;
using namespace exageostat::configurations::data_configurations;
using namespace std;

/**
 * @brief The main function of the program.
 *
 * This function generates synthetic data for ExaGeoStat using the provided command line arguments.
 *
 * @param argc The number of command line arguments.
 * @param argv The command line arguments.
 * @return The status code of the program.
 */

int main(int argc, char **argv) {
    // Create a unique pointer to a DataGenerator object
    unique_ptr<DataGenerator> synthetic_generator;

    // Create a new synthetic_data_configurations object with the provided command line arguments
    SyntheticDataConfigurations* synthetic_data_configurations = SyntheticDataConfigurations::GetInstance();
    synthetic_data_configurations->InitializeArguments(argc, argv);

    // Create the DataGenerator object
    synthetic_generator = synthetic_generator->CreateGenerator(synthetic_data_configurations);

    // Initialize the locations of the generated data
    synthetic_generator->GenerateLocations();

    // Define a struct to hold pointers to the x, y, and z coordinates of the generated data
    struct DataPointers {
        double *x;
        double *y;
        double *z;
    } dataPointers{};

    // Set the pointers in the DataPointers struct to the location coordinates of the generated data
    dataPointers.x = synthetic_generator->GetLocations()->GetLocationX();
    dataPointers.y = synthetic_generator->GetLocations()->GetLocationY();
    dataPointers.z = synthetic_generator->GetLocations()->GetLocationZ();

    // Print the generated location coordinates
    cout << "Generated Locations are .. " << endl;
    int timeSlot;
    if (synthetic_data_configurations->GetDimension() != DimensionST) {
        timeSlot = 1;
    } else {
        timeSlot = synthetic_data_configurations->GetTimeSlot();
    }
    for (auto i = 0; i < synthetic_data_configurations->GetProblemSize() * timeSlot; i++) {
        cout << "X: " << dataPointers.x[i] << " Y: " << dataPointers.y[i];
        if (synthetic_data_configurations->GetDimension() != Dimension2D) {
            cout << " Z: " << dataPointers.z[i];
        }
        cout << endl;
    }

    synthetic_generator->GenerateDescriptors();
    synthetic_generator->GenerateObservations();

    return 0;
}