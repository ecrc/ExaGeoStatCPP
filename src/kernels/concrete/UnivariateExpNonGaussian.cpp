
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file UnivariateExpNonGaussian.cpp
 *
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-04-14
**/

#include <kernels/concrete/UnivariateExpNonGaussian.hpp>

using namespace exageostat::kernels;
using namespace exageostat::dataunits;
using namespace std;

UnivariateExpNonGaussian::UnivariateExpNonGaussian() {
    this->mP = 1;
    this->mParametersNumber = 6;
}

Kernel *UnivariateExpNonGaussian::Create() {
    return new UnivariateExpNonGaussian();
}

namespace exageostat::kernels {
    bool UnivariateExpNonGaussian::plugin_name = plugins::PluginRegistry<exageostat::kernels::Kernel>::Add(
            "UnivariateExpNonGaussian", UnivariateExpNonGaussian::Create);
}

void UnivariateExpNonGaussian::GenerateCovarianceMatrix(double *apMatrixA, int aRowsNumber, int aColumnsNumber,
                                                            int aRowOffset, int aColumnOffset, Locations *apLocation1,
                                                            Locations *apLocation2, Locations *apLocation3,
                                                            double *apLocalTheta, int aDistanceMetric) {

    int i, j;
    int i0 = aRowOffset;
    int j0 = aColumnOffset;
    double x0, y0, z0;
    double expr = 0.0;
    double sigma_square = 1;

    for (i = 0; i < aRowsNumber; i++) {
        j0 = aColumnOffset;
        for (j = 0; j < aColumnsNumber; j++) {
            expr = CalculateDistance(apLocation1, apLocation2, i0, j0, aDistanceMetric, 0) / apLocalTheta[0];

            if (expr == 0)
                apMatrixA[i + j * aRowsNumber] = sigma_square /*+ 1e-4*/;
            else
                apMatrixA[i + j * aRowsNumber] = exp(-expr);

            j0++;
        }
        i0++;
    }
}