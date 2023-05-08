
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file UnivariateMaternStationary.cpp
 *
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-04-14
**/

#include <kernels/concrete/UnivariateMaternStationary.hpp>

using namespace exageostat::kernels;
using namespace exageostat::dataunits;
using namespace std;

UnivariateMaternStationary::UnivariateMaternStationary() {
    this->mP = 1;
}

Kernel *UnivariateMaternStationary::Create() {
    return new UnivariateMaternStationary();
}

void UnivariateMaternStationary::GenerateCovarianceMatrix(double *apMatrixA, int aRowsNumber, int aColumnsNumber,
                                                          int aRowOffset, int aColumnOffset, Locations *apLocation1,
                                                          Locations *apLocation2, Locations *apLocation3,
                                                          double *apLocalTheta, int aDistanceMetric) {
    const double sigma_square = apLocalTheta[0];
    const double nu = apLocalTheta[2];
    const double inv_con = 1.0 / (sigma_square * pow(2, nu - 1) * tgamma(nu));
    int i0 = aRowOffset;

    for (int i = 0; i < aRowsNumber; i++) {
        int j0 = aColumnOffset;
        for (int j = 0; j < aColumnsNumber; j++) {
            const double dist =
                    CalculateDistance(apLocation1, apLocation2, j0, i0, aDistanceMetric, 0) / apLocalTheta[1];
            *(apMatrixA + i + j * aRowsNumber) = (dist == 0.0)
                                                 ? sigma_square
                                                 : inv_con * pow(dist, nu) * gsl_sf_bessel_Knu(nu, dist);
            j0++;
        }
        i0++;
    }
}

namespace exageostat::kernels {
    bool UnivariateMaternStationary::plugin_name = plugins::PluginRegistry<exageostat::kernels::Kernel>::Add(
            "UnivariateMaternStationary", UnivariateMaternStationary::Create);
}
