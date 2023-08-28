
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file UnivariateMaternNonStationary.cpp
 * @brief Implementation of the UnivariateMaternNonStationary kernel.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @author Mahmoud ElKarargy
 * @date 2023-04-13
**/

#include <kernels/concrete/UnivariateMaternNonStationary.hpp>

using namespace std;

using namespace exageostat::kernels;
using namespace exageostat::dataunits;

template<typename T>
UnivariateMaternNonStationary<T>::UnivariateMaternNonStationary() {
    this->mP = 1;
    this->mParametersNumber = 9;
}

template<typename T>
Kernel<T> *UnivariateMaternNonStationary<T>::Create() {
    return new UnivariateMaternNonStationary();
}

namespace exageostat::kernels {
    template<typename T> bool UnivariateMaternNonStationary<T>::plugin_name = plugins::PluginRegistry<exageostat::kernels::Kernel<T>>::Add(
            "UnivariateMaternNonStationary", UnivariateMaternNonStationary<T>::Create);
}

template<typename T>
void UnivariateMaternNonStationary<T>::GenerateCovarianceMatrix(T *apMatrixA, const int &aRowsNumber,
                                                                const int &aColumnsNumber,
                                                                const int &aRowOffset, const int &aColumnOffset,
                                                                dataunits::Locations<T> &aLocation1,
                                                                dataunits::Locations<T> &aLocation2,
                                                                dataunits::Locations<T> &aLocation3, T *aLocalTheta,
                                                                const int &aDistanceMetric) {

    double location1X, location1Y, location2X, location2Y, location3X, location3Y;
    double theta_0i, theta_0j, theta_1i, theta_1j, theta_2i, theta_2j;
    double dx, dy;
    double dist;
    double con, sigma_square, beta, nu;
    int i, j;

    aLocation3 = aLocation1;
    double x_max = aLocation1->GetLocationX()[0];
    double x_min = aLocation1->GetLocationX()[0];
    double y_max = aLocation1->GetLocationY()[0];
    double y_min = aLocation1->GetLocationY()[0];
    for (i = 1; i < 9; i++) {
        if (x_max < aLocation1->GetLocationX()[i])
            x_max = aLocation1->GetLocationX()[i];
        if (x_min > aLocation1->GetLocationX()[i])
            x_min = aLocation1->GetLocationX()[i];
        if (y_max < aLocation1->GetLocationY()[i])
            y_max = aLocation1->GetLocationY()[i];
        if (y_max > aLocation1->GetLocationY()[i])
            y_max = aLocation1->GetLocationY()[i];
    }

    aLocation3->GetLocationX()[0] = x_min + (x_max - x_min) / 2;
    aLocation3->GetLocationY()[0] = y_min + (y_max - y_min) / 2;
    printf(" The central point is ( %f, %f)\n", aLocation3->GetLocationX()[0], aLocation3->GetLocationY()[0]);

    // Compute the covariance matrix elements
    for (j = 0; j < aColumnsNumber; j++) {
        location1X = aLocation1->GetLocationX()[aColumnOffset + j];
        location1Y = aLocation1->GetLocationY()[aColumnOffset + j];
        location3X = aLocation3->GetLocationX()[j];
        location3Y = aLocation3->GetLocationY()[j];
        dx = abs(location1X - location3X);
        dy = abs(location1Y - location3Y);

        theta_0i = aLocalTheta[0] + (aLocalTheta[1] * dx) + (aLocalTheta[2] * dy);
        theta_1i = aLocalTheta[3] + (aLocalTheta[4] * dx) + (aLocalTheta[5] * dy);
        theta_2i = aLocalTheta[6] + (aLocalTheta[7] * dx) + (aLocalTheta[8] * dy);

        for (i = 0; i < aRowsNumber; i++) {
            location2X = aLocation2->GetLocationX()[aRowOffset + i];
            location2Y = aLocation2->GetLocationY()[aRowOffset + i];
            location3X = aLocation3->GetLocationX()[i];
            location3Y = aLocation3->GetLocationY()[i];
            dx = abs(location2X - location3X);
            dy = abs(location2Y - location3Y);

            theta_0j = aLocalTheta[0] + (aLocalTheta[1] * dx) + (aLocalTheta[2] * dy);
            theta_1j = aLocalTheta[3] + (aLocalTheta[4] * dx) + (aLocalTheta[5] * dy);
            theta_2j = aLocalTheta[6] + (aLocalTheta[7] * dx) + (aLocalTheta[8] * dy);

            // Compute the Matérn parameters and distance metric
            sigma_square = (theta_0i + theta_0j) / 2;
            beta = (theta_1i + theta_1j) / 2;
            nu = (theta_2i + theta_2j) / 2;
            con = pow(2, (nu - 1)) / tgamma(nu);

            con = 1.0 / con;
            con = sigma_square * con;
            int flag = 0;

            //MLE calculation
            dist = CalculateDistance(aLocation1, aLocation2, i, j, aDistanceMetric, flag) / beta;

            *(apMatrixA + i + j * aRowsNumber) = (dist == 0.0)
                                                 ? sigma_square
                                                 //                                                 : inv_con * pow(dist, nu) * gsl_sf_bessel_Knu(nu, dist);
                                                 : con * pow(dist, nu) * gsl_sf_bessel_Knu(nu, dist);
        }
    }
}
