
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file UnivariateSpacetimeMaternStationary.cpp
 *
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-04-14
**/
#include <kernels/concrete/UnivariateSpacetimeMaternStationary.hpp>
#include<cmath>
#include <gsl/gsl_sf_bessel.h>

using namespace exageostat::kernels;
using namespace exageostat::dataunits;
using namespace std;


void UnivariateSpacetimeMaternStationary::GenerateCovarianceMatrix(double *apMatrixA, int aRowsNumber, int aColumnsNumber,
                                                                 int aRowOffset, int aColumnOffset, Locations *apLocation1,
                                                                 Locations *apLocation2, Locations *apLocation3,
                                                                 double *apLocalTheta, int aDistanceMetric) {
    int i, j;
    int i0 = aRowOffset;
    int j0 = aColumnOffset;
    double x0, y0, z0, z1;
    double expr = 0.0, expr1 = 0.0, expr2 = 0.0, expr3 = 0.0, expr4 = 0.0;
    double con = 0.0;
    double sigma_square = apLocalTheta[0];

    con = pow(2, (apLocalTheta[2] - 1)) * tgamma(apLocalTheta[2]);
    con = 1.0 / con;
    con = sigma_square * con;

    for (i = 0; i < aRowsNumber; i++) {
        j0 = aColumnOffset;
        z0 = apLocation1->GetLocationZ()[i0];
        for (j = 0; j < aColumnsNumber; j++) {
            z1 = apLocation2->GetLocationZ()[j0];

            expr = CalculateDistance(apLocation1, apLocation2, i0, j0, aDistanceMetric, 1) / apLocalTheta[1];
            expr2 = pow(pow(sqrt(pow(z0 - z1, 2)), 2 * apLocalTheta[4]) / apLocalTheta[3] + 1.0, apLocalTheta[5] / 2.0);
            expr3 = expr / expr2;
            expr4 = pow(pow(sqrt(pow(z0 - z1, 2)), 2 * apLocalTheta[4]) / apLocalTheta[3] + 1.0,
                        apLocalTheta[5] + apLocalTheta[6]);

            if (expr == 0)
                apMatrixA[i + j * aRowsNumber] = sigma_square / expr4 /*+ 1e-4*/;
            else
                apMatrixA[i + j * aRowsNumber] = con * pow(expr3, apLocalTheta[2])
                                                 * gsl_sf_bessel_Knu(apLocalTheta[2], expr3) / expr4; // Matern Function
            j0++;
        }
        i0++;
    }
}