
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file BivariateSpacetimeMaternStationary.cpp
 *
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-04-14
**/

#include <kernels/concrete/BivariateSpacetimeMaternStationary.hpp>

using namespace exageostat::kernels;
using namespace exageostat::dataunits;
using namespace std;

BivariateSpacetimeMaternStationary::BivariateSpacetimeMaternStationary() {
    this->mP = 2;
    this->mParametersNumber = 10;
}

Kernel *BivariateSpacetimeMaternStationary::Create() {
    return new BivariateSpacetimeMaternStationary();
}

namespace exageostat::kernels {
    bool BivariateSpacetimeMaternStationary::plugin_name = plugins::PluginRegistry<exageostat::kernels::Kernel>::Add(
            "BivariateSpacetimeMaternStationary", BivariateSpacetimeMaternStationary::Create);
}

void BivariateSpacetimeMaternStationary::GenerateCovarianceMatrix(double *apMatrixA, int aRowsNumber, int aColumnsNumber,
                                                               int aRowOffset, int aColumnOffset, Locations *apLocation1,
                                                               Locations *apLocation2, Locations *apLocation3,
                                                               double *apLocalTheta, int aDistanceMetric) {
    int i, j;
    int i0 = aRowOffset;
    int j0 = aColumnOffset;
    double x0, y0, z0, z1;
    double expr = 0.0, expr1 = 0.0, expr2 = 0.0, expr3 = 0.0, expr4 = 0.0;
    double con1 = 0.0, con2 = 0.0, con12 = 0.0, rho = 0.0, nu12 = 0.0;

    con1 = pow(2, (apLocalTheta[3] - 1)) * tgamma(apLocalTheta[3]);
    con1 = 1.0 / con1;
    con1 = apLocalTheta[0] * con1;

    con2 = pow(2, (apLocalTheta[4] - 1)) * tgamma(apLocalTheta[4]);
    con2 = 1.0 / con2;
    con2 = apLocalTheta[1] * con2;

    //The average
    nu12 = 0.5 * (apLocalTheta[3] + apLocalTheta[4]);

    rho = apLocalTheta[5] * sqrt((tgamma(apLocalTheta[3] + 1) * tgamma(apLocalTheta[4] + 1)) /
                                 (tgamma(apLocalTheta[3]) * tgamma(apLocalTheta[4]))) *
          tgamma(nu12) / tgamma(nu12 + 1);

    con12 = pow(2, (nu12 - 1)) * tgamma(nu12);
    con12 = 1.0 / con12;
    con12 = rho * sqrt(apLocalTheta[0] * apLocalTheta[1]) * con12;

    i0 /= 2;
    for (i = 0; i < aRowsNumber - 1; i += 2) {
        j0 = aColumnOffset / 2;
        z0 = apLocation1->GetLocationZ()[i0];

        for (j = 0; j < aColumnsNumber - 1; j += 2) {
            z1 = apLocation2->GetLocationZ()[j0];

            expr = CalculateDistance(apLocation1, apLocation2, i0, j0, aDistanceMetric, 1) / (apLocalTheta[2] * 1000);
            expr2 = pow(pow(sqrt(pow(z0 - z1, 2)), 2 * apLocalTheta[7]) / apLocalTheta[6] + 1, apLocalTheta[8] / 2);
            expr3 = expr / expr2;
            expr4 = pow(pow(sqrt(pow(z0 - z1, 2)), 2 * apLocalTheta[7]) / apLocalTheta[6] + 1,
                        apLocalTheta[8] + apLocalTheta[9]);

            if (expr == 0) {
                apMatrixA[i + j * aRowsNumber] = apLocalTheta[0] / expr4;
                apMatrixA[(i + 1) + j * aRowsNumber] = apMatrixA[i + (j + 1) * aRowsNumber] = rho
                                                                                              * sqrt(apLocalTheta[0] * apLocalTheta[1]) / expr4;
                apMatrixA[(i + 1) + (j + 1) * aRowsNumber] = apLocalTheta[1] / expr4;
            } else {
                apMatrixA[i + j * aRowsNumber] = con1 * pow(expr3, apLocalTheta[3])
                                                 * gsl_sf_bessel_Knu(apLocalTheta[3], expr3) / expr4;
                apMatrixA[(i + 1) + j * aRowsNumber] = apMatrixA[i + (j + 1) * aRowsNumber] = con12 * pow(expr3, nu12)
                                                                                              * gsl_sf_bessel_Knu(nu12, expr3) / expr4;
                apMatrixA[(i + 1) + (j + 1) * aRowsNumber] = con2 * pow(expr3, apLocalTheta[4])
                                                             * gsl_sf_bessel_Knu(apLocalTheta[4], expr3) / expr4;
            }
            j0++;
        }
        i0++;
    }
}