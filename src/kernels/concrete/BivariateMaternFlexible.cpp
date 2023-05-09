
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file BivariateMaternFlexible.cpp
 *
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-04-14
**/

#include <kernels/concrete/BivariateMaternFlexible.hpp>

using namespace exageostat::kernels;
using namespace exageostat::dataunits;
using namespace std;

BivariateMaternFlexible::BivariateMaternFlexible() {
    this->mP = 1;
}

Kernel *BivariateMaternFlexible::Create() {
    return new BivariateMaternFlexible();
}

void BivariateMaternFlexible::GenerateCovarianceMatrix(double *apMatrixA, int aRowsNumber, int aColumnsNumber,
                                                          int aRowOffset, int aColumnOffset, Locations *apLocation1,
                                                          Locations *apLocation2, Locations *apLocation3,
                                                          double *apLocalTheta, int aDistanceMetric) {
    int i, j;
    int i0 = aRowOffset;
    int j0 = aColumnOffset;
    double x0, y0;
    double expr1 = 0.0, expr2 = 0.0, expr12 = 0.0;
    double con1 = 0.0, con2 = 0.0, con12 = 0.0, scale12 = 0.0, rho = 0.0, nu12 = 0.0, sigma_square11 = 0.0, sigma_square22 = 0.0;
    double scale1 = apLocalTheta[0], scale2 = apLocalTheta[1], nu1 = apLocalTheta[4], nu2 = apLocalTheta[5];

    scale12 = pow(0.5 * (pow(scale1, -2) + pow(scale2, -2)) + apLocalTheta[2] * (1 - apLocalTheta[3]),
                  -0.5); //Remark 1 (c) of Apanasovich et al. (2012)

    nu12 = 0.5 * (nu1 + nu2) + apLocalTheta[6] * (1 - apLocalTheta[7]); //Theorem 1 (i) of Apanasovich et al. (2012).

    rho = apLocalTheta[8] * apLocalTheta[9] * apLocalTheta[10] *
          pow(scale12, 2 * apLocalTheta[6] + (nu1 + nu2))
          * tgamma(0.5 * (nu1 + nu2) + 1) * tgamma(nu12) /
          tgamma(nu12 + 1); //Equation (8) of Apanasovich et al. (2012).

    sigma_square11 = apLocalTheta[8] * apLocalTheta[8] *
                     pow(scale1, 2 * apLocalTheta[6] + nu1 + nu1) *
                     tgamma(nu1); //Equation (8) of Apanasovich et al. (2012).

    sigma_square22 = apLocalTheta[9] * apLocalTheta[9] *
                     pow(scale2, 2 * apLocalTheta[6] + nu2 + nu2) *
                     tgamma(nu2); //Equation (8) of Apanasovich et al. (2012).

    con1 = pow(2, (nu1 - 1)) * tgamma(nu1);
    con1 = 1.0 / con1;
    con1 = sigma_square11 * con1;

    con2 = pow(2, (nu2 - 1)) * tgamma(nu2);
    con2 = 1.0 / con2;
    con2 = sigma_square22 * con2;

    con12 = pow(2, (nu12 - 1)) * tgamma(nu12);
    con12 = 1.0 / con12;
    con12 = rho * con12;

    i0 /= 2;
    for (i = 0; i < aRowsNumber; i += 2) {
        j0 = aColumnOffset / 2;
        for (j = 0; j < aColumnsNumber; j += 2) {
            expr1 = CalculateDistance(apLocation1, apLocation2, i0, j0, aDistanceMetric, 0) / scale1;
            expr2 = CalculateDistance(apLocation1, apLocation2, i0, j0, aDistanceMetric, 0) / scale2;
            expr12 = CalculateDistance(apLocation1, apLocation2, i0, j0, aDistanceMetric, 0) / scale12;

            if (expr1 == 0) {
                apMatrixA[i + j * aRowsNumber] = apLocalTheta[0];
                apMatrixA[(i + 1) + j * aRowsNumber] = apMatrixA[i + (j + 1) * aRowsNumber] = rho;
                apMatrixA[(i + 1) + (j + 1) * aRowsNumber] = apLocalTheta[1];
            } else {
                apMatrixA[i + j * aRowsNumber] = con1 * pow(expr1, nu1)
                                                 * gsl_sf_bessel_Knu(nu1, expr1);
                apMatrixA[(i + 1) + j * aRowsNumber] = apMatrixA[i + (j + 1) * aRowsNumber] = con12
                                                                                              * pow(expr12, nu12) * gsl_sf_bessel_Knu(nu12, expr12);
                apMatrixA[(i + 1) + (j + 1) * aRowsNumber] = con2 * pow(expr2, nu2)
                                                             * gsl_sf_bessel_Knu(nu2, expr2);
            }
            j0++;
        }
        i0++;
    }
}

namespace exageostat::kernels {
    bool BivariateMaternFlexible::plugin_name = plugins::PluginRegistry<exageostat::kernels::Kernel>::Add(
            "BivariateMaternFlexible", BivariateMaternFlexible::Create);
}