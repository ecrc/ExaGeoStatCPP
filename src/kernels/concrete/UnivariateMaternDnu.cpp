
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file UnivariateMaternDnu.cpp
 * @brief Implementation of the UnivariateMaternDnu kernel.
 * @version 1.0.1
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-04-14
**/

#include <kernels/concrete/UnivariateMaternDnu.hpp>


using namespace exageostat::kernels;
using namespace exageostat::dataunits;
using namespace exageostat::helpers;

template<typename T>
UnivariateMaternDnu<T>::UnivariateMaternDnu() {
    this->mP = 1;
    this->mParametersNumber = 3;
}

template<typename T>
Kernel<T> *UnivariateMaternDnu<T>::Create() {
    KernelsConfigurations::GetParametersNumberKernelMap()["UnivariateMaternDnu"] = 3;
    return new UnivariateMaternDnu();
}

namespace exageostat::kernels {
    template<typename T> bool UnivariateMaternDnu<T>::plugin_name = plugins::PluginRegistry<exageostat::kernels::Kernel<T>>::Add(
            "UnivariateMaternDnu", UnivariateMaternDnu<T>::Create);
}

template<typename T>
void UnivariateMaternDnu<T>::GenerateCovarianceMatrix(T *apMatrixA, const int &aRowsNumber, const int &aColumnsNumber,
                                                      const int &aRowOffset, const int &aColumnOffset,
                                                      Locations<T> &aLocation1, Locations<T> &aLocation2,
                                                      Locations<T> &aLocation3, T *aLocalTheta,
                                                      const int &aDistanceMetric) {

    int i, j;
    int i0 = aRowOffset;
    int j0;
    T expr;
    T nu_expr;
    T sigma_square = aLocalTheta[0];
    int flag = aLocation1.GetLocationZ() == nullptr ? 0 : 1;

    for (i = 0; i < aRowsNumber; i++) {
        j0 = aColumnOffset;
        for (j = 0; j < aColumnsNumber; j++) {
            expr = DistanceCalculationHelpers<T>::CalculateDistance(aLocation1, aLocation2, i0, j0, aDistanceMetric,
                                                                    flag) / aLocalTheta[1];
            if (expr == 0) {
                apMatrixA[i + j * aRowsNumber] = 0.0;
            } else {
                //derivative with respect to nu
                nu_expr = -2 * log(2.0) * pow(2, -aLocalTheta[2]) * 1 / tgamma(aLocalTheta[2]) *
                          pow(expr, aLocalTheta[2]) * gsl_sf_bessel_Knu(aLocalTheta[2], expr) +
                          pow(2, 1 - aLocalTheta[2]) *
                          (-1 / tgamma(aLocalTheta[2]) * gsl_sf_psi(aLocalTheta[2]) * pow(expr, aLocalTheta[2]) *
                           gsl_sf_bessel_Knu(aLocalTheta[2], expr) + 1 / tgamma(aLocalTheta[2]) *
                                                                     (pow(expr, aLocalTheta[2]) * log(expr) *
                                                                      gsl_sf_bessel_Knu(aLocalTheta[2], expr) +
                                                                      pow(expr, aLocalTheta[2]) *
                                                                      BasselFunction<T>::CalculateDerivativeBesselNu(aLocalTheta[2], expr)));
                apMatrixA[i + j * aRowsNumber] = sigma_square * nu_expr;
            }
            j0++;

        }
        i0++;
    }
}