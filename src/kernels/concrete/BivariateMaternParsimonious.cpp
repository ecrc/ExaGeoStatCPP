
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file BivariateMaternParsimonious.cpp
 * @brief Implementation of the BivariateMaternParsimonious kernel.
 * @version 1.0.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-04-14
**/

#include <kernels/concrete/BivariateMaternParsimonious.hpp>
#include <helpers/DistanceCalculationHelpers.hpp>

using namespace exageostat::kernels;
using namespace exageostat::dataunits;
using namespace exageostat::helpers;

template<typename T>
BivariateMaternParsimonious<T>::BivariateMaternParsimonious() {
    this->mP = 2;
    this->mParametersNumber = 6;
}

template<typename T>
Kernel<T> *BivariateMaternParsimonious<T>::Create() {
    KernelsConfigurations::GetParametersNumberKernelMap()["BivariateMaternParsimonious"] = 6;
    return new BivariateMaternParsimonious();
}

namespace exageostat::kernels {
    template<typename T> bool BivariateMaternParsimonious<T>::plugin_name = plugins::PluginRegistry<exageostat::kernels::Kernel<T>>::Add(
            "BivariateMaternParsimonious", BivariateMaternParsimonious::Create);
}

template<typename T>
void BivariateMaternParsimonious<T>::GenerateCovarianceMatrix(T *apMatrixA, const int &aRowsNumber,
                                                              const int &aColumnsNumber, const int &aRowOffset,
                                                              const int &aColumnOffset, Locations<T> &aLocation1,
                                                              Locations<T> &aLocation2, Locations<T> &aLocation3,
                                                              T *aLocalTheta, const int &aDistanceMetric) {

    int i, j;
    int i0 = aRowOffset;
    int j0;
    T expr;
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
    int flag = aLocation1.GetLocationZ() == nullptr ? 0 : 1;

    for (i = 0; i < aRowsNumber; i += 2) {
        j0 = aColumnOffset / 2;
        for (j = 0; j < aColumnsNumber; j += 2) {
            expr = DistanceCalculationHelpers<T>::CalculateDistance(aLocation1, aLocation2, i0, j0, aDistanceMetric,
                                                                    0) / aLocalTheta[2];
            if (expr == 0) {
                apMatrixA[i + j * aRowsNumber] = aLocalTheta[0];

                if (((i + 1) + j * aRowsNumber) < aRowsNumber * aColumnsNumber) {
                    apMatrixA[(i + 1) + j * aRowsNumber] = rho * sqrt(aLocalTheta[0] * aLocalTheta[1]);
                }
                if ((i + (j + 1) * aRowsNumber) < aRowsNumber * aColumnsNumber) {
                    apMatrixA[i + (j + 1) * aRowsNumber] = rho * sqrt(aLocalTheta[0] * aLocalTheta[1]);
                }
                if (((i + 1) + (j + 1) * aRowsNumber) < aRowsNumber * aColumnsNumber) {
                    apMatrixA[(i + 1) + (j + 1) * aRowsNumber] = aLocalTheta[1];
                }

            } else {
                apMatrixA[i + j * aRowsNumber] =
                        con1 * pow(expr, aLocalTheta[3]) * gsl_sf_bessel_Knu(aLocalTheta[3], expr);

                if (((i + 1) + j * aRowsNumber) < aRowsNumber * aColumnsNumber) {
                    apMatrixA[(i + 1) + j * aRowsNumber] = con12 * pow(expr, nu12) * gsl_sf_bessel_Knu(nu12, expr);
                }
                if ((i + (j + 1) * aRowsNumber) < aRowsNumber * aColumnsNumber) {
                    apMatrixA[i + (j + 1) * aRowsNumber] = con12 * pow(expr, nu12) * gsl_sf_bessel_Knu(nu12, expr);
                }
                if (((i + 1) + (j + 1) * aRowsNumber) < aRowsNumber * aColumnsNumber) {
                    apMatrixA[(i + 1) + (j + 1) * aRowsNumber] =
                            con2 * pow(expr, aLocalTheta[4]) * gsl_sf_bessel_Knu(aLocalTheta[4], expr);
                }
            }
            j0++;
        }
        i0++;
    }
}
