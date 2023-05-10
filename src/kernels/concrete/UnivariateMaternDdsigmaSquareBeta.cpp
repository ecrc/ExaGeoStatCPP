
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file UnivariateMaternDdsigmaSquareBeta.cpp
 *
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-04-14
**/

#include <kernels/concrete/UnivariateMaternDdsigmaSquareBeta.hpp>

using namespace exageostat::kernels;
using namespace exageostat::dataunits;
using namespace std;

UnivariateMaternDdsigmaSquareBeta::UnivariateMaternDdsigmaSquareBeta() {
    /// TODO: FIX THEIR VALUES
    this->mP = 1;
    this->mParametersNumber = 3;
}

Kernel *UnivariateMaternDdsigmaSquareBeta::Create() {
    return new UnivariateMaternDdsigmaSquareBeta();
}

namespace exageostat::kernels {
    bool UnivariateMaternDdsigmaSquareBeta::plugin_name = plugins::PluginRegistry<exageostat::kernels::Kernel>::Add(
            "UnivariateMaternDdsigmaSquareBeta", UnivariateMaternDdsigmaSquareBeta::Create);
}

void UnivariateMaternDdsigmaSquareBeta::GenerateCovarianceMatrix(double *apMatrixA, int aRowsNumber, int aColumnsNumber,
                                                                 int aRowOffset, int aColumnOffset,
                                                                 Locations *apLocation1,
                                                                 Locations *apLocation2, Locations *apLocation3,
                                                                 double *apLocalTheta, int aDistanceMetric) {
    int i, j;
    int i0 = aRowOffset;
    int j0 = aColumnOffset;
    double x0, y0, z0;
    double expr = 0.0;
    double con = 0.0;
    double beta_expr = 0.0;
    double sigma_square = apLocalTheta[0];
    con = pow(2, (apLocalTheta[2] - 1)) * tgamma(apLocalTheta[2]);
    con = 1.0 / con;
    for (i = 0; i < aRowsNumber; i++) {
        j0 = aColumnOffset;
        for (j = 0; j < aColumnsNumber; j++) {
            expr = CalculateDistance(apLocation1, apLocation2, i0, j0, aDistanceMetric, 0) / apLocalTheta[1];
            if (expr == 0) {

                apMatrixA[i + j * aRowsNumber] = 0.0;

            } else {
                beta_expr = -apLocalTheta[2] / apLocalTheta[1] * pow(expr, apLocalTheta[2])
                            * gsl_sf_bessel_Knu(apLocalTheta[2], expr) - pow(expr, apLocalTheta[2])
                                                                         * (apLocalTheta[2] / expr *
                                                                            gsl_sf_bessel_Knu(apLocalTheta[2], expr) -
                                                                            gsl_sf_bessel_Knu(apLocalTheta[2] + 1,
                                                                                              expr)) * expr /
                                                                         apLocalTheta[1];
                apMatrixA[i + j * aRowsNumber] = con * beta_expr; // derivative with respect to sigma square and beta
            }
            j0++;
        }
        i0++;
    }
}