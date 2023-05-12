
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file BivariateMaternParsimonious.cpp
 *
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-04-14
**/

#include <kernels/concrete/BivariateMaternParsimonious.hpp>
#include<cmath>
#include <gsl/gsl_sf_bessel.h>

using namespace exageostat::kernels;
using namespace exageostat::dataunits;
using namespace std;

BivariateMaternParsimonious::BivariateMaternParsimonious() {
    this->mP = 2;
    this->mParametersNumber = 6;
}

Kernel *BivariateMaternParsimonious::Create() {
    return new BivariateMaternParsimonious();
}

namespace exageostat::kernels {
    bool BivariateMaternParsimonious::plugin_name = plugins::PluginRegistry<exageostat::kernels::Kernel>::Add(
            "BivariateMaternParsimonious", BivariateMaternParsimonious::Create);
}

void BivariateMaternParsimonious::GenerateCovarianceMatrix(double *apMatrixA, int aRowsNumber, int aColumnsNumber,
                                                          int aRowOffset, int aColumnOffset, Locations *apLocation1,
                                                          Locations *apLocation2, Locations *apLocation3,
                                                          double *apLocalTheta, int aDistanceMetric) {
    int i, j;
    int i0 = aRowOffset;
    int j0 = aColumnOffset;
    double expr = 0.0;
    double con1, con2, con12, rho, nu12;

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

    printf("con1: %f con2: %f nu12: %f rho: %f con12: %f\n", con1, con2, nu12, rho, con12);

    i0 /= 2;
    for (i = 0; i < aRowsNumber - 1; i += 2) {
        j0 = aColumnOffset / 2;
        for (j = 0; j < aColumnsNumber - 1; j += 2) {
            expr = CalculateDistance(apLocation1, apLocation2, i0, j0, aDistanceMetric, 0) / apLocalTheta[2];

            if (expr == 0) {
                apMatrixA[i + j * aRowsNumber] = apLocalTheta[0];
                apMatrixA[(i + 1) + j * aRowsNumber] = apMatrixA[i + (j + 1) * aRowsNumber] = rho
                                                                                              * sqrt(apLocalTheta[0] * apLocalTheta[1]);
                apMatrixA[(i + 1) + (j + 1) * aRowsNumber] = apLocalTheta[1];
            } else {
                apMatrixA[i + j * aRowsNumber] = con1 * pow(expr, apLocalTheta[3])
                                                 * gsl_sf_bessel_Knu(apLocalTheta[3], expr);
                apMatrixA[(i + 1) + j * aRowsNumber] = apMatrixA[i + (j + 1) * aRowsNumber] = con12 * pow(expr, nu12)
                                                                                              * gsl_sf_bessel_Knu(nu12, expr);
                apMatrixA[(i + 1) + (j + 1) * aRowsNumber] = con2 * pow(expr, apLocalTheta[4])
                                                             * gsl_sf_bessel_Knu(apLocalTheta[4], expr);
            }
            j0++;
        }
        i0++;
    }
    std::cout << "after: " << std::endl;

    for (j = 0; j < aColumnsNumber; j++) {
        for (i = 0; i < aRowsNumber; i++) {
            std::cout << *(apMatrixA + i + j * aRowsNumber) << " ";
        }
        std::cout << std::endl;
    }
}
