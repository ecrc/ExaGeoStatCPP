
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file UnivariateMaternDdsigmaSquareNu.cpp
 *
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-04-14
**/

#include <kernels/concrete/UnivariateMaternDdsigmaSquareNu.hpp>

using namespace exageostat::kernels;
using namespace exageostat::dataunits;
using namespace std;

UnivariateMaternDdsigmaSquareNu::UnivariateMaternDdsigmaSquareNu() {
    /// TODO: FIX THEIR VALUES
    this->mP = 1;
    this->mParametersNumber = 3;
}

Kernel *UnivariateMaternDdsigmaSquareNu::Create() {
    return new UnivariateMaternDdsigmaSquareNu();
}

namespace exageostat::kernels {
    bool UnivariateMaternDdsigmaSquareNu::plugin_name = plugins::PluginRegistry<exageostat::kernels::Kernel>::Add(
            "UnivariateMaternDdsigmaSquareNu", UnivariateMaternDdsigmaSquareNu::Create);
}

void UnivariateMaternDdsigmaSquareNu::GenerateCovarianceMatrix(double *apMatrixA, int aRowsNumber, int aColumnsNumber,
                                                   int aRowOffset, int aColumnOffset, Locations *apLocation1,
                                                   Locations *apLocation2, Locations *apLocation3,
                                                   double *apLocalTheta, int aDistanceMetric) {
    int i, j;
    int i0 = aRowOffset;
    int j0 = aColumnOffset;
    double x0, y0, z0;
    double expr = 0.0;
    double con = 0.0;
    double nu_expr = 0.0;
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
                //// TODO: Make sure gsl_sf_bessel_Kn have the same effect as psi in C
                nu_expr = (1 - apLocalTheta[2]) * 1 / pow(2, apLocalTheta[2]) * 1 / tgamma(apLocalTheta[2])
                          * pow(expr, apLocalTheta[2]) * gsl_sf_bessel_Knu(apLocalTheta[2], expr) +
                          pow(2, 1 - apLocalTheta[2])
                          * (-1 / tgamma(apLocalTheta[2]) * gsl_sf_psi(apLocalTheta[2]) * pow(expr, apLocalTheta[2])
                             * gsl_sf_bessel_Knu(apLocalTheta[2], expr) + 1 / tgamma(apLocalTheta[2])
                                                                          * (pow(expr, apLocalTheta[2]) * log(expr)
                                                                             * gsl_sf_bessel_Knu(apLocalTheta[2], expr) +
                                                                             pow(expr, apLocalTheta[2])
                                                                             * (gsl_sf_bessel_Knu(apLocalTheta[2] + 0.000000001, expr) - gsl_sf_bessel_Knu(apLocalTheta[2], expr)) / 0.000000001));
                apMatrixA[i + j * aRowsNumber] = nu_expr; // derivative with respect to sigma square
            }
            j0++;
        }
        i0++;
    }
}