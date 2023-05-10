
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file TrivariateMaternParsimonious.cpp
 *
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-04-14
**/

#include <kernels/concrete/TrivariateMaternParsimonious.hpp>

using namespace exageostat::kernels;
using namespace exageostat::dataunits;
using namespace std;

TrivariateMaternParsimonious::TrivariateMaternParsimonious() {
    this->mP = 3;
    this->mParametersNumber = 10;
}

Kernel *TrivariateMaternParsimonious::Create() {
    return new TrivariateMaternParsimonious();
}

namespace exageostat::kernels {
    bool TrivariateMaternParsimonious::plugin_name = plugins::PluginRegistry<exageostat::kernels::Kernel>::Add(
            "TrivariateMaternParsimonious", TrivariateMaternParsimonious::Create);
}

void TrivariateMaternParsimonious::GenerateCovarianceMatrix(double *apMatrixA, int aRowsNumber, int aColumnsNumber,
                                                                  int aRowOffset, int aColumnOffset, Locations *apLocation1,
                                                                  Locations *apLocation2, Locations *apLocation3,
                                                                  double *apLocalTheta, int aDistanceMetric) {
    int i, j;
    int i0 = aRowOffset;
    int j0 = aColumnOffset;
    double x0, y0;
    double expr = 0.0;
    double con1 = 0.0, con2 = 0.0, con3 = 0.0, con12 = 0.0,
            con13 = 0.0, con23 = 0.0, rho12 = 0.0, rho13 = 0.0, rho23 = 0.0, nu12 = 0.0, nu13 = 0.0, nu23 = 0.0;

    con1 = pow(2, (apLocalTheta[4] - 1)) * tgamma(apLocalTheta[4]);
    con1 = 1.0 / con1;
    con1 = apLocalTheta[0] * con1;

    con2 = pow(2, (apLocalTheta[5] - 1)) * tgamma(apLocalTheta[5]);
    con2 = 1.0 / con2;
    con2 = apLocalTheta[1] * con2;

    con3 = pow(2, (apLocalTheta[6] - 1)) * tgamma(apLocalTheta[6]);
    con3 = 1.0 / con3;
    con3 = apLocalTheta[2] * con3;

    // The average
    nu12 = 0.5 * (apLocalTheta[4] + apLocalTheta[5]);
    nu13 = 0.5 * (apLocalTheta[4] + apLocalTheta[6]);
    nu23 = 0.5 * (apLocalTheta[5] + apLocalTheta[6]);

    rho12 = apLocalTheta[7] * sqrt((tgamma(apLocalTheta[4] + 1) * tgamma(apLocalTheta[5] + 1)) /
                                   (tgamma(apLocalTheta[4]) * tgamma(apLocalTheta[5]))) *
            tgamma(nu12) / tgamma(nu12 + 1);

    rho13 = apLocalTheta[8] * sqrt((tgamma(apLocalTheta[4] + 1) * tgamma(apLocalTheta[6] + 1)) /
                                   (tgamma(apLocalTheta[4]) * tgamma(apLocalTheta[6]))) *
            tgamma(nu13) / tgamma(nu13 + 1);

    rho23 = apLocalTheta[9] * sqrt((tgamma(apLocalTheta[5] + 1) * tgamma(apLocalTheta[6] + 1)) /
                                   (tgamma(apLocalTheta[5]) * tgamma(apLocalTheta[6]))) *
            tgamma(nu23) / tgamma(nu23 + 1);

    con12 = pow(2, (nu12 - 1)) * tgamma(nu12);
    con12 = 1.0 / con12;
    con12 = rho12 * sqrt(apLocalTheta[0] * apLocalTheta[1]) * con12;

    con13 = pow(2, (nu13 - 1)) * tgamma(nu13);
    con13 = 1.0 / con13;
    con13 = rho13 * sqrt(apLocalTheta[0] * apLocalTheta[2]) * con13;

    con23 = pow(2, (nu23 - 1)) * tgamma(nu23);
    con23 = 1.0 / con23;
    con23 = rho23 * sqrt(apLocalTheta[1] * apLocalTheta[2]) * con23;

    i0 /= 3;
    for (i = 0; i < aRowsNumber - 1; i += 3) {
        j0 = aColumnOffset / 3;
        for (j = 0; j < aColumnsNumber - 1; j += 3) {
            expr = CalculateDistance(apLocation1, apLocation2, i0, j0, aDistanceMetric, 0) / apLocalTheta[3];

            if (expr == 0) {
                apMatrixA[i + j * aRowsNumber] = apLocalTheta[0];

                apMatrixA[(i + 1) + j * aRowsNumber] = apMatrixA[i + (j + 1) * aRowsNumber] = rho12 * sqrt(apLocalTheta[0] * apLocalTheta[1]);
                apMatrixA[(i + 2) + j * aRowsNumber] = apMatrixA[i + (j + 2) * aRowsNumber] = rho13 * sqrt(apLocalTheta[0] * apLocalTheta[2]);

                apMatrixA[(i + 1) + (j + 1) * aRowsNumber] = apLocalTheta[1];

                apMatrixA[(i + 1) + (j + 2) * aRowsNumber] = apMatrixA[(i + 2) + (j + 1) * aRowsNumber] = rho23 * sqrt(apLocalTheta[1] * apLocalTheta[2]);

                apMatrixA[(i + 2) + (j + 2) * aRowsNumber] = apLocalTheta[2];
            } else {
                apMatrixA[i + j * aRowsNumber] = con1 * pow(expr, apLocalTheta[4]) * gsl_sf_bessel_Knu(apLocalTheta[4], expr);

                apMatrixA[(i + 1) + j * aRowsNumber] = apMatrixA[i + (j + 1) * aRowsNumber] = con12 * pow(expr, nu12) * gsl_sf_bessel_Knu(nu12, expr);
                apMatrixA[(i + 2) + j * aRowsNumber] = apMatrixA[i + (j + 2) * aRowsNumber] = con13 * pow(expr, nu13) * gsl_sf_bessel_Knu(nu13, expr);

                apMatrixA[(i + 1) + (j + 1) * aRowsNumber] = con2 * pow(expr, apLocalTheta[5]) * gsl_sf_bessel_Knu(apLocalTheta[5], expr);

                apMatrixA[(i + 1) + (j + 2) * aRowsNumber] = apMatrixA[(i + 2) + (j + 1) * aRowsNumber] = con23 * pow(expr, nu23) * gsl_sf_bessel_Knu(nu23, expr);

                apMatrixA[(i + 2) + (j + 2) * aRowsNumber] = con3 * pow(expr, apLocalTheta[6]) * gsl_sf_bessel_Knu(apLocalTheta[6], expr);
            }
            j0++;
        }
        i0++;
    }
}
