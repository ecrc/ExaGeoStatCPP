
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file BivariateSpacetimeMaternStationary.cpp
 * @brief Implementation of the BivariateSpacetimeMaternStationary kernel.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-04-14
**/

#include <kernels/concrete/BivariateSpacetimeMaternStationary.hpp>
#include <helpers/DistanceCalculationHelpers.hpp>

using namespace exageostat::kernels;
using namespace exageostat::dataunits;
using namespace exageostat::helpers;

template<typename T>
BivariateSpacetimeMaternStationary<T>::BivariateSpacetimeMaternStationary() {
    this->mP = 2;
    this->mParametersNumber = 10;
}

template<typename T>
Kernel<T> *BivariateSpacetimeMaternStationary<T>::Create() {
    KernelsConfigurations::GetParametersNumberKernelMap()["BivariateSpacetimeMaternStationary"] = 10;
    return new BivariateSpacetimeMaternStationary();
}

namespace exageostat::kernels {
    template<typename T> bool BivariateSpacetimeMaternStationary<T>::plugin_name = plugins::PluginRegistry<exageostat::kernels::Kernel<T>>::Add(
            "BivariateSpacetimeMaternStationary", BivariateSpacetimeMaternStationary<T>::Create);
}

template<typename T>
void
BivariateSpacetimeMaternStationary<T>::GenerateCovarianceMatrix(T *apMatrixA, const int &aRowsNumber,
                                                                const int &aColumnsNumber, const int &aRowOffset,
                                                                const int &aColumnOffset, Locations<T> &aLocation1,
                                                                Locations<T> &aLocation2, Locations<T> &aLocation3,
                                                                T *aLocalTheta, const int &aDistanceMetric) {

    int i, j;
    int i0 = aRowOffset;
    int j0;
    T z0, z1;
    T expr, expr2, expr3, expr4;
    T con1, con2, con12, rho, nu12;

    con1 = pow(2, (aLocalTheta[3] - 1)) * tgamma(aLocalTheta[3]);
    con1 = 1.0 / con1;
    con1 = aLocalTheta[0] * con1;
    con2 = pow(2, (aLocalTheta[4] - 1)) * tgamma(aLocalTheta[4]);
    con2 = 1.0 / con2;
    con2 = aLocalTheta[1] * con2;

    //The average
    nu12 = 0.5 * (aLocalTheta[3] + aLocalTheta[4]);
    rho = aLocalTheta[5] * sqrt((tgamma(aLocalTheta[3] + 1) * tgamma(aLocalTheta[4] + 1)) /
                                (tgamma(aLocalTheta[3]) * tgamma(aLocalTheta[4]))) * tgamma(nu12) / tgamma(nu12 + 1);
    con12 = pow(2, (nu12 - 1)) * tgamma(nu12);
    con12 = 1.0 / con12;
    con12 = rho * sqrt(aLocalTheta[0] * aLocalTheta[1]) * con12;

    i0 /= 2;
    int matrix_size = aRowsNumber * aColumnsNumber;
    int index;
    int flag = 1;

    for (i = 0; i < aRowsNumber; i += 2) {
        j0 = aColumnOffset / 2;
        if (aLocation1.GetLocationZ() != nullptr) {
            z0 = aLocation1.GetLocationZ()[i0];
        }

        for (j = 0; j < aColumnsNumber; j += 2) {
            if (aLocation2.GetLocationZ() != nullptr) {
                z1 = aLocation2.GetLocationZ()[j0];
            }
            expr = DistanceCalculationHelpers<T>::CalculateDistance(aLocation1, aLocation2, i0, j0, aDistanceMetric,
                                                                    flag) / (aLocalTheta[2] * 1000);
            expr2 = pow(pow(sqrt(pow(z0 - z1, 2)), 2 * aLocalTheta[7]) / aLocalTheta[6] + 1, aLocalTheta[8] / 2);
            expr3 = expr / expr2;
            expr4 = pow(pow(sqrt(pow(z0 - z1, 2)), 2 * aLocalTheta[7]) / aLocalTheta[6] + 1,
                        aLocalTheta[8] + aLocalTheta[9]);

            if (expr == 0) {
                apMatrixA[i + j * aRowsNumber] = aLocalTheta[0] / expr4;
                index = (i + 1) + j * aRowsNumber;
                if (index < matrix_size) {
                    apMatrixA[(i + 1) + j * aRowsNumber] = rho * sqrt(aLocalTheta[0] * aLocalTheta[1]) / expr4;
                }
                index = i + (j + 1) * aRowsNumber;
                if (index < matrix_size) {
                    apMatrixA[i + (j + 1) * aRowsNumber] = rho * sqrt(aLocalTheta[0] * aLocalTheta[1]) / expr4;
                }
                index = (i + 1) + (j + 1) * aRowsNumber;
                if (index < matrix_size) {
                    apMatrixA[(i + 1) + (j + 1) * aRowsNumber] = aLocalTheta[1] / expr4;
                }
            } else {
                apMatrixA[i + j * aRowsNumber] =
                        con1 * pow(expr3, aLocalTheta[3]) * gsl_sf_bessel_Knu(aLocalTheta[3], expr3) / expr4;
                index = (i + 1) + j * aRowsNumber;
                if (index < matrix_size) {
                    apMatrixA[(i + 1) + j * aRowsNumber] =
                            con12 * pow(expr3, nu12) * gsl_sf_bessel_Knu(nu12, expr3) / expr4;
                }
                index = i + (j + 1) * aRowsNumber;
                if (index < matrix_size) {
                    apMatrixA[i + (j + 1) * aRowsNumber] =
                            con12 * pow(expr3, nu12) * gsl_sf_bessel_Knu(nu12, expr3) / expr4;
                }
                index = (i + 1) + (j + 1) * aRowsNumber;
                if (index < matrix_size) {
                    apMatrixA[(i + 1) + (j + 1) * aRowsNumber] =
                            con2 * pow(expr3, aLocalTheta[4]) * gsl_sf_bessel_Knu(aLocalTheta[4], expr3) / expr4;
                }
            }
            j0++;
        }
        i0++;
    }
}