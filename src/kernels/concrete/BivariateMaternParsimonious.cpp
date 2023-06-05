
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
                                                           double *aLocalTheta, int aDistanceMetric) {
    int i, j;
    int i0 = aRowOffset;
    int j0;
    double expr;
    double con1, con2, con12, rho, nu12;

    con1 = pow(2, (aLocalTheta[3] - 1)) * tgamma(aLocalTheta[3]);
    con1 = 1.0 / con1;
    con1 = aLocalTheta[0] * con1;

    con2 = pow(2, (aLocalTheta[4] - 1)) * tgamma(aLocalTheta[4]);
    con2 = 1.0 / con2;
    con2 = aLocalTheta[1] * con2;

    //The average
    nu12 = 0.5 * (aLocalTheta[3] + aLocalTheta[4]);

    rho = aLocalTheta[5] * sqrt((tgamma(aLocalTheta[3] + 1) * tgamma(aLocalTheta[4] + 1)) /
                                (tgamma(aLocalTheta[3]) * tgamma(aLocalTheta[4]))) *
          tgamma(nu12) / tgamma(nu12 + 1);


    con12 = pow(2, (nu12 - 1)) * tgamma(nu12);
    con12 = 1.0 / con12;
    con12 = rho * sqrt(aLocalTheta[0] * aLocalTheta[1]) * con12;

    i0 /= 2;
    for (i = 0; i < aRowsNumber; i += 2) {
        j0 = aColumnOffset / 2;
        for (j = 0; j < aColumnsNumber; j += 2) {
            expr = CalculateDistance(apLocation1, apLocation2, i0, j0, aDistanceMetric, 0) / aLocalTheta[2];
            if (expr == 0) {
                apMatrixA[i + j * aRowsNumber] = aLocalTheta[0];

                if(((i + 1) + j * aRowsNumber ) < aRowsNumber * aColumnsNumber){
                    apMatrixA[(i + 1) + j * aRowsNumber] = rho * sqrt(aLocalTheta[0] * aLocalTheta[1]);
                }
                if((i + (j + 1) * aRowsNumber) < aRowsNumber * aColumnsNumber){
                    apMatrixA[i + (j + 1) * aRowsNumber] = rho * sqrt(aLocalTheta[0] * aLocalTheta[1]);
                }
                if(((i + 1) + (j + 1) * aRowsNumber) < aRowsNumber * aColumnsNumber){
                    apMatrixA[(i + 1) + (j + 1) * aRowsNumber] = aLocalTheta[1];
                }

            } else {
                apMatrixA[i + j * aRowsNumber] = con1 * pow(expr, aLocalTheta[3])
                                                 * gsl_sf_bessel_Knu(aLocalTheta[3], expr);

                if(((i + 1) + j * aRowsNumber ) < aRowsNumber * aColumnsNumber){
                    apMatrixA[(i + 1) + j * aRowsNumber] = con12 * pow(expr, nu12) * gsl_sf_bessel_Knu(nu12,expr);
                }
                if((i + (j + 1) * aRowsNumber) < aRowsNumber * aColumnsNumber){
                    apMatrixA[i + (j + 1) * aRowsNumber] = con12 * pow(expr, nu12) * gsl_sf_bessel_Knu(nu12,expr);
                }
                if(((i + 1) + (j + 1) * aRowsNumber) < aRowsNumber * aColumnsNumber){
                    apMatrixA[(i + 1) + (j + 1) * aRowsNumber] = con2 * pow(expr, aLocalTheta[4]) * gsl_sf_bessel_Knu(aLocalTheta[4], expr);
                }
            }
            j0++;
        }
        i0++;
    }

}
