
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file UnivariateSpacetimeMaternStationary.cpp
 * @brief Implementation of the UnivariateSpacetimeMaternStationary kernel.
 * @version 1.0.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-04-14
**/

#include<cmath>

#include <gsl/gsl_sf_bessel.h>

#include <kernels/concrete/UnivariateSpacetimeMaternStationary.hpp>

using namespace exageostat::kernels;
using namespace exageostat::dataunits;

template<typename T>
UnivariateSpacetimeMaternStationary<T>::UnivariateSpacetimeMaternStationary() {
    this->mP = 1;
    this->mParametersNumber = 7;
}

template<typename T>
Kernel<T> *UnivariateSpacetimeMaternStationary<T>::Create() {
    KernelsConfigurations::GetParametersNumberKernelMap()["UnivariateSpacetimeMaternStationary"] = 7;
    return new UnivariateSpacetimeMaternStationary();
}

namespace exageostat::kernels {
    template<typename T> bool UnivariateSpacetimeMaternStationary<T>::plugin_name = plugins::PluginRegistry<exageostat::kernels::Kernel<T>>::Add(
            "UnivariateSpacetimeMaternStationary", UnivariateSpacetimeMaternStationary<T>::Create);
}

template<typename T>
void
UnivariateSpacetimeMaternStationary<T>::GenerateCovarianceMatrix(T *apMatrixA, const int &aRowsNumber,
                                                                 const int &aColumnsNumber, const int &aRowOffset,
                                                                 const int &aColumnOffset, Locations<T> &aLocation1,
                                                                 Locations<T> &aLocation2, Locations<T> &aLocation3,
                                                                 T *aLocalTheta, const int &aDistanceMetric) {

    int i, j;
    int i0 = aRowOffset;
    int j0;
    double z0, z1;
    double expr, expr2, expr3, expr4;
    double con;
    double sigma_square = aLocalTheta[0];

    con = pow(2, (aLocalTheta[2] - 1)) * tgamma(aLocalTheta[2]);
    con = 1.0 / con;
    con = sigma_square * con;
    int flag = 1;

    for (i = 0; i < aRowsNumber; i++) {
        j0 = aColumnOffset;
        z0 = aLocation1.GetLocationZ()[i0];
        for (j = 0; j < aColumnsNumber; j++) {
            z1 = aLocation2.GetLocationZ()[j0];

            expr = this->CalculateDistance(aLocation1, aLocation2, i0, j0, aDistanceMetric, flag) / aLocalTheta[1];
            expr2 = pow(pow(sqrt(pow(z0 - z1, 2)), 2 * aLocalTheta[4]) / aLocalTheta[3] + 1.0, aLocalTheta[5] / 2.0);
            expr3 = expr / expr2;
            expr4 = pow(pow(sqrt(pow(z0 - z1, 2)), 2 * aLocalTheta[4]) / aLocalTheta[3] + 1.0,
                        aLocalTheta[5] + aLocalTheta[6]);

            if (expr == 0) {
                apMatrixA[i + j * aRowsNumber] = sigma_square / expr4;
            } else {
                // Matern Function
                apMatrixA[i + j * aRowsNumber] =
                        con * pow(expr3, aLocalTheta[2]) * gsl_sf_bessel_Knu(aLocalTheta[2], expr3) / expr4;
            }
            j0++;
        }
        i0++;
    }
}