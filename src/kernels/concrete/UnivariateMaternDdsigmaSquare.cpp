
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file UnivariateMaternDdsigmaSquare.cpp
 *
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-04-14
**/

#include <kernels/concrete/UnivariateMaternDdsigmaSquare.hpp>

using namespace exageostat::kernels;
using namespace exageostat::dataunits;
using namespace std;

UnivariateMaternDdsigmaSquare::UnivariateMaternDdsigmaSquare() {
    /// TODO: FIX THEIR VALUES
    this->mP = 1;
    this->mParametersNumber = 3;
}

Kernel *UnivariateMaternDdsigmaSquare::Create() {
    return new UnivariateMaternDdsigmaSquare();
}

namespace exageostat::kernels {
    bool UnivariateMaternDdsigmaSquare::plugin_name = plugins::PluginRegistry<exageostat::kernels::Kernel>::Add(
            "UnivariateMaternDdsigmaSquare", UnivariateMaternDdsigmaSquare::Create);
}

void UnivariateMaternDdsigmaSquare::GenerateCovarianceMatrix(double *apMatrixA, int aRowsNumber, int aColumnsNumber,
                                                             int aRowOffset, int aColumnOffset, Locations *apLocation1,
                                                             Locations *apLocation2, Locations *apLocation3,
                                                             double *apLocalTheta, int aDistanceMetric) {
    int i, j;
    //// TODO: Implementation is Empty in the old version!
    for (i = 0; i < aRowsNumber; i++) {
        for (j = 0; j < aColumnsNumber; j++) {
            apMatrixA[i + j * aRowsNumber] = 0.0;
        }
    }
}