/*
 * Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
 * Copyright (C) 2023 by Brightskies inc,
 * All rights reserved.
 * ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).
 */

/**
* @file synthatic_data.cpp
* @version 1.0.0
* @author Sameh Abdulah
* @date 2023-01-31
**/

#include <iostream>
#include <configurations/data-generation/concrete/SyntheticDataConfigurations.hpp>

using namespace std;
using namespace exageostat::configurations::data_configurations;

int main(int argc, char **argv) {

    auto *syntheticDataConfigurations = new SyntheticDataConfigurations(argc, argv);

    int N = syntheticDataConfigurations->GetProblemSize();
    if (N != 0) {
        cout << "You set N by: " << N << endl;
    }

    string kernel = syntheticDataConfigurations->GetKernel();
    if (!kernel.empty()) {
        cout << "You set Kernel by: " << kernel << endl;
    }

    string dimension = syntheticDataConfigurations->GetDimension();
    if (!dimension.empty()) {
        cout << "You set Dimension by: " << dimension << endl;
    }

    int PGrid = syntheticDataConfigurations->GetPGrid();
    if (PGrid != 0) {
        cout << "You set N by: " << PGrid << endl;
    }

    delete syntheticDataConfigurations;
    return 0;
}
