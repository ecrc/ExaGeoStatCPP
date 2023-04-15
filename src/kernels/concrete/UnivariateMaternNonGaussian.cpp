
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file UnivariateMaternNonGaussian.cpp
 *
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-04-14
**/

#include <kernels/concrete/UnivariateMaternNonGaussian.hpp>
#include<cmath>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_psi.h>


using namespace exageostat::kernels;
using namespace exageostat::dataunits;
using namespace std;


void UnivariateMaternNonGaussian::GenerateCovarianceMatrix(double *apMatrixA, int aRowsNumber, int aColumnsNumber,
                                                      int aRowOffset, int aColumnOffset, Locations *apLocation1,
                                                      Locations *apLocation2, Locations *apLocation3,
                                                      double *apLocalTheta, int aDistanceMetric) {
    //localtheta[0] <- \phi
    //localtheta[1] <- \nu
    int i, j;
    int i0 = aRowOffset;
    int j0 = aColumnOffset;
    double x0, y0, z0;
    double expr = 0.0;
    double con = 0.0;
    double sigma_square = 1;

    con = pow(2, (apLocalTheta[1] - 1)) * tgamma(apLocalTheta[1]);
    con = 1.0 / con;
    con = sigma_square * con;

    for (i = 0; i < aRowsNumber; i++) {
        j0 = aColumnOffset;
        for (j = 0; j < aColumnsNumber; j++) {
            expr = 4 * sqrt(2 * apLocalTheta[1]) *
                   (CalculateDistance(apLocation1, apLocation2, i0, j0, aDistanceMetric, 0) / apLocalTheta[0]);
            if (expr == 0)
                apMatrixA[i + j * aRowsNumber] = sigma_square /*+ 1e-4*/;
            else
                apMatrixA[i + j * aRowsNumber] = con * pow(expr, apLocalTheta[1])
                                                 * gsl_sf_bessel_Knu(apLocalTheta[1], expr); // Matern Function
            j0++;
        }
        i0++;
    }
}