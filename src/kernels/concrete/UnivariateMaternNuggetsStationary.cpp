
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file UnivariateMaternNuggetsStationary.cpp
 *
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-04-14
**/

#include <kernels/concrete/UnivariateMaternNuggetsStationary.hpp>
#include<cmath>
#include <gsl/gsl_sf_bessel.h>

using namespace exageostat::kernels;
using namespace exageostat::dataunits;
using namespace std;


void UnivariateMaternNuggetsStationary::GenerateCovarianceMatrix(double *apMatrixA, int aRowsNumber, int aColumnsNumber,
                                                       int aRowOffset, int aColumnOffset, Locations *apLocation1,
                                                       Locations *apLocation2, Locations *apLocation3,
                                                       double *apLocalTheta, int aDistanceMetric) {
    int i, j;
    int i0 = aRowOffset;
    int j0 = aColumnOffset;
    double x0, y0, z0;
    double expr = 0.0;
    double con = 0.0;
    double sigma_square = apLocalTheta[0];// * apLocalTheta[0];

    con = pow(2, (apLocalTheta[2] - 1)) * tgamma(apLocalTheta[2]);
    con = 1.0 / con;
    con = sigma_square * con;
    if (apLocation1->GetLocationZ() == nullptr || apLocation2->GetLocationZ() == nullptr) {
        for (i = 0; i < aRowsNumber; i++) {
            j0 = aColumnOffset;
            for (j = 0; j < aColumnsNumber; j++) {
                expr = CalculateDistance(apLocation1, apLocation2, j0, i0, aDistanceMetric, 0) / apLocalTheta[1];
                if (expr == 0)
                    apMatrixA[i + j * aRowsNumber] = sigma_square + apLocalTheta[3];
                else
                    apMatrixA[i + j * aRowsNumber] =
                            con * pow(expr, apLocalTheta[2]) * gsl_sf_bessel_Knu(apLocalTheta[2], expr); // Matern Function
                j0++;
            }
            i0++;
        }
    } else {
        for (i = 0; i < aRowsNumber; i++) {
            j0 = aColumnOffset;
            for (j = 0; j < aColumnsNumber; j++) {
                expr = CalculateDistance(apLocation1, apLocation2, j0, i0, aDistanceMetric, 1);
                if (expr == 0)
                    apMatrixA[i + j * aRowsNumber] = sigma_square + apLocalTheta[3];
                else
                    apMatrixA[i + j * aRowsNumber] =
                            con * pow(expr, apLocalTheta[2]) * gsl_sf_bessel_Knu(apLocalTheta[2], expr); // Matern Function
                j0++;
            }
            i0++;
        }
    }
}