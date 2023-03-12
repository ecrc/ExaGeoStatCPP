
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file SyntheticGenerator.cpp
 *
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-03-04
**/

#include <iostream>
#include <data-generators/DataGenerator.hpp>

using namespace exageostat::generators;
using namespace exageostat::dataunits;
using namespace std;
using namespace exageostat::configurations::data_configurations;

int main(int argc, char **argv) {
    unique_ptr<DataGenerator> syntheticGenerator;

    // Object has automatic storage duration (usually is on the stack)
    auto syntheticDataConfigurations = new SyntheticDataConfigurations(argc, argv);

    syntheticGenerator = syntheticGenerator->CreateGenerator(syntheticDataConfigurations);
    syntheticGenerator->InitializeLocations();

    struct DataPointers{
        double *x;
        double *y;
        double *z;
    }dataPointers;

    dataPointers.x = syntheticGenerator->GetLocations()->GetLocationX();
    dataPointers.y = syntheticGenerator->GetLocations()->GetLocationY();
    dataPointers.z = syntheticGenerator->GetLocations()->GetLocationZ();

    cout << "Generated Locations are .. " << endl;
    int timeSlot;
    if (syntheticDataConfigurations->GetDimension() != DimensionST){
        timeSlot = 1;
    }
    else {
        timeSlot = syntheticDataConfigurations->GetTimeSlot();
    }
    for(auto i = 0; i < syntheticDataConfigurations->GetProblemSize() * timeSlot; i ++){
        cout << "X: " << dataPointers.x[i] << " Y: " << dataPointers.y[i];
        if (syntheticDataConfigurations->GetDimension() != Dimension2D){
            cout << " Z: " << dataPointers.z[i];
        }
        cout << endl;
    }

    return 0;
}