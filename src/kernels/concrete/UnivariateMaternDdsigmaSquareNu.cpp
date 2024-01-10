
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file UnivariateMaternDdsigmaSquareNu.cpp
 * @brief Implementation of the UnivariateMaternDdsigmaSquareNu kernel.
 * @version 1.0.1
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-04-14
**/

#include <kernels/concrete/UnivariateMaternDdsigmaSquareNu.hpp>
#include <helpers/DistanceCalculationHelpers.hpp>

using namespace exageostat::kernels;
using namespace exageostat::dataunits;
using namespace exageostat::helpers;

template<typename T>
UnivariateMaternDdsigmaSquareNu<T>::UnivariateMaternDdsigmaSquareNu() {
    this->mP = 1;
    this->mParametersNumber = 3;
}

template<typename T>
Kernel<T> *UnivariateMaternDdsigmaSquareNu<T>::Create() {
    KernelsConfigurations::GetParametersNumberKernelMap()["UnivariateMaternDdsigmaSquareNu"] = 3;
    return new UnivariateMaternDdsigmaSquareNu();
}

namespace exageostat::kernels {
    template<typename T> bool UnivariateMaternDdsigmaSquareNu<T>::plugin_name = plugins::PluginRegistry<exageostat::kernels::Kernel<T>>::Add(
            "UnivariateMaternDdsigmaSquareNu", UnivariateMaternDdsigmaSquareNu<T>::Create);
}

template<typename T>
void UnivariateMaternDdsigmaSquareNu<T>::GenerateCovarianceMatrix(T *apMatrixA, const int &aRowsNumber,
                                                                  const int &aColumnsNumber, const int &aRowOffset,
                                                                  const int &aColumnOffset, Locations<T> &aLocation1,
                                                                  Locations<T> &aLocation2, Locations<T> &aLocation3,
                                                                  T *aLocalTheta, const int &aDistanceMetric) {

    int i, j;
    int i0 = aRowOffset;
    int j0;
    T expr;
    T nu_expr;
    int flag = aLocation1.GetLocationZ() == nullptr ? 0 : 1;

    for (i = 0; i < aRowsNumber; i++) {
        j0 = aColumnOffset;
        for (j = 0; j < aColumnsNumber; j++) {
            expr = DistanceCalculationHelpers<T>::CalculateDistance(aLocation1, aLocation2, i0, j0, aDistanceMetric,
                                                                    flag) / aLocalTheta[1];
            if (expr == 0) {
                apMatrixA[i + j * aRowsNumber] = 0.0;
            } else {
                // derivative with respect to sigma square
                nu_expr = (1 - aLocalTheta[2]) * 1 / pow(2, aLocalTheta[2]) * 1 / tgamma(aLocalTheta[2]) *
                          pow(expr, aLocalTheta[2]) * gsl_sf_bessel_Knu(aLocalTheta[2], expr) +
                          pow(2, 1 - aLocalTheta[2]) *
                          (-1 / tgamma(aLocalTheta[2]) * gsl_sf_psi(aLocalTheta[2]) * pow(expr, aLocalTheta[2]) *
                           gsl_sf_bessel_Knu(aLocalTheta[2], expr) + 1 / tgamma(aLocalTheta[2]) *
                                                                     (pow(expr, aLocalTheta[2]) * log(expr) *
                                                                      gsl_sf_bessel_Knu(aLocalTheta[2], expr) +
                                                                      pow(expr, aLocalTheta[2]) *
                                                                      (gsl_sf_bessel_Knu(aLocalTheta[2] + 0.000000001,
                                                                                         expr) -
                                                                       gsl_sf_bessel_Knu(aLocalTheta[2], expr)) /
                                                                      0.000000001));
                apMatrixA[i + j * aRowsNumber] = nu_expr;
            }
            j0++;
        }
        i0++;
    }
}