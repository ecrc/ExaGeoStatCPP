
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file UnivariateMaternNonStationary.cpp
 *
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-04-13
**/

#include <kernels/concrete/UnivariateMaternNonStationary.hpp>
#include <iostream>
#include<cmath>
#include <gsl/gsl_sf_bessel.h>

using namespace exageostat::kernels;
using namespace exageostat::dataunits;
using namespace std;


void UnivariateMaternStationary::GenerateCovarianceMatrix(double *apMatrixA, int aRowsNumber, int aColumnsNumber,
                                                          int aRowOffset, int aColumnOffset, Locations *apLocation1,
                                                          Locations *apLocation2, Locations *apLocation3,
                                                          double *apLocalTheta, int aDistanceMetric) {
    double location1X, location1Y, location2X, location2Y, location3X, location3Y;
    double theta_0i, theta_0j, theta_1i, theta_1j, theta_2i, theta_2j;
    double dx, dy;
    double expr = 0.0;
    double con, sigma_square, beta, nu;
    int i, j;

    // Compute the covariance matrix elements
    for (j = 0; j < aColumnsNumber; j++) {
        location1X = apLocation1->GetLocationX()[aColumnOffset + j];
        location1Y = apLocation1->GetLocationY()[aColumnOffset + j];
        location3X = apLocation3->GetLocationX()[j];
        location3Y = apLocation3->GetLocationY()[j];

        dx = abs(location1X - location3X);
        dy = abs(location1Y - location3Y);
        theta_0i = apLocalTheta[0] + (apLocalTheta[1] * dx) + (apLocalTheta[2] * dy);
        theta_1i = apLocalTheta[3] + (apLocalTheta[4] * dx) + (apLocalTheta[5] * dy);
        theta_2i = apLocalTheta[6] + (apLocalTheta[7] * dx) + (apLocalTheta[8] * dy);

        for (i = 0; i < aRowsNumber; i++) {
            location2X = apLocation2->GetLocationX()[aRowOffset + i];
            location2Y = apLocation2->GetLocationY()[aRowOffset + i];
            location3X = apLocation3->GetLocationX()[i];
            location3Y = apLocation3->GetLocationY()[i];

            dx = abs(location2X - location3X);
            dy = abs(location2Y - location3Y);
            theta_0j = apLocalTheta[0] + (apLocalTheta[1] * dx) + (apLocalTheta[2] * dy);
            theta_1j = apLocalTheta[3] + (apLocalTheta[4] * dx) + (apLocalTheta[5] * dy);
            theta_2j = apLocalTheta[6] + (apLocalTheta[7] * dx) + (apLocalTheta[8] * dy);

            // Compute the Matérn parameters and distance metric
            sigma_square = (theta_0i + theta_0j) / 2;
            beta = (theta_1i + theta_1j) / 2;
            nu = (theta_2i + theta_2j) / 2;
            con = pow(2, (nu - 1)) / tgamma(nu);

            con = 1.0 / con;
            con = sigma_square * con;
            //MLE calculation
            expr = CalculateDistance(apLocation1, apLocation2, i, j, aDistanceMetric, 0) / beta;
            if (expr == 0)
                apMatrixA[i + j * aRowsNumber] = sigma_square; /* + 1e-4*/
            else
                apMatrixA[i + j * aRowsNumber] = con * pow(expr, nu) * gsl_sf_bessel_Knu(nu, expr); // Matern Function

//// TODO: What do you think?
//            dx = abs(location1X - location2X);
//            dy = abs(location1Y - location2Y);
//
//            // Compute the covariance matrix element based on the distance metric
//            switch (aDistanceMetric) {
//                case 1: // Euclidean distance
//                    expr = beta * dx;
//                    break;
//                case 2: // Manhattan distance
//                    expr = beta * sqrt(pow(dx, 2) + pow(dy, 2));
//                    break;
//                case 3: // Minkowski distance
//                    expr = beta * pow(pow(dx, nu) + pow(dy, nu), (1 / nu));
//                    expr *= con;
//                    break;
//            }
//            apMatrixA[i + j * aRowsNumber] = sigma_square * exp(-expr);
        }
    }
}