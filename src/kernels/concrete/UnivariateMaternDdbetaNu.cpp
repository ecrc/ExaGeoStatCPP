
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file UnivariateMaternDdbetaNu.cpp
 * @brief Implementation of the UnivariateMaternDdbetaNu kernel.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-04-14
**/

#include <kernels/concrete/UnivariateMaternDdbetaNu.hpp>
#include <helpers/DistanceCalculationHelpers.hpp>

using namespace exageostat::kernels;
using namespace exageostat::dataunits;
using namespace exageostat::helpers;

template<typename T>
UnivariateMaternDdbetaNu<T>::UnivariateMaternDdbetaNu() {
    this->mP = 1;
    this->mParametersNumber = 3;
}

template<typename T>
Kernel<T> *UnivariateMaternDdbetaNu<T>::Create() {
    KernelsConfigurations::GetParametersNumberKernelMap()["UnivariateMaternDdbetaNu"] = 3;
    return new UnivariateMaternDdbetaNu();
}

namespace exageostat::kernels {
    template<typename T> bool UnivariateMaternDdbetaNu<T>::plugin_name = plugins::PluginRegistry<exageostat::kernels::Kernel<T>>::Add(
            "UnivariateMaternDdbetaNu", UnivariateMaternDdbetaNu<T>::Create);
}

template<typename T>
void
UnivariateMaternDdbetaNu<T>::GenerateCovarianceMatrix(T *apMatrixA, const int &aRowsNumber, const int &aColumnsNumber,
                                                      const int &aRowOffset, const int &aColumnOffset,
                                                      Locations<T> &aLocation1, Locations<T> &aLocation2,
                                                      Locations<T> &aLocation3, T *aLocalTheta,
                                                      const int &aDistanceMetric) {

    int i, j;
    int i0 = aRowOffset;
    int j0;
    T expr;
    T con;
    T nu_expr;
    T nu_expr_prime;
    T sigma_square = aLocalTheta[0];
    con = pow(2, (aLocalTheta[2] - 1)) * tgamma(aLocalTheta[2]);
    con = 1.0 / con;
    int flag = aLocation1.GetLocationZ() == nullptr ? 0 : 1;

    for (i = 0; i < aRowsNumber; i++) {
        j0 = aColumnOffset;
        for (j = 0; j < aColumnsNumber; j++) {
            expr = DistanceCalculationHelpers<T>::CalculateDistance(aLocation1, aLocation2, i0, j0, aDistanceMetric,
                                                                    flag) / aLocalTheta[1];
            if (expr == 0) {
                apMatrixA[i + j * aRowsNumber] = 0.0;
            } else {
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
                nu_expr_prime = (1 - aLocalTheta[2]) * 1 / pow(2, aLocalTheta[2]) * 1 / tgamma(aLocalTheta[2]) *
                                pow(expr, aLocalTheta[2]) *
                                (aLocalTheta[2] / expr * gsl_sf_bessel_Knu(aLocalTheta[2], expr) -
                                 gsl_sf_bessel_Knu(aLocalTheta[2] + 1, expr)) + pow(2, 1 - aLocalTheta[2]) *
                                                                                (-1 / tgamma(aLocalTheta[2]) *
                                                                                 gsl_sf_psi(aLocalTheta[2]) *
                                                                                 pow(expr, aLocalTheta[2]) *
                                                                                 (aLocalTheta[2] / expr *
                                                                                  gsl_sf_bessel_Knu(aLocalTheta[2],
                                                                                                    expr) -
                                                                                  gsl_sf_bessel_Knu(aLocalTheta[2] + 1,
                                                                                                    expr)) +
                                                                                 1 / tgamma(aLocalTheta[2]) *
                                                                                 (pow(expr, aLocalTheta[2]) *
                                                                                  log(expr) * (aLocalTheta[2] / expr *
                                                                                               gsl_sf_bessel_Knu(
                                                                                                       aLocalTheta[2],
                                                                                                       expr) -
                                                                                               gsl_sf_bessel_Knu(
                                                                                                       aLocalTheta[2] +
                                                                                                       1, expr)) +
                                                                                  pow(expr, aLocalTheta[2]) *
                                                                                  this->CalculateSecondDerivativeBesselNuInput(
                                                                                          aLocalTheta[2], expr)));
                apMatrixA[i + j * aRowsNumber] = (-1 / aLocalTheta[1] * (con * pow(expr, aLocalTheta[2]) *
                                                                         gsl_sf_bessel_Knu(aLocalTheta[2], expr)) -
                                                  aLocalTheta[2] / aLocalTheta[1] * nu_expr -
                                                  expr / aLocalTheta[1] * nu_expr_prime) * sigma_square;
            }
            j0++;
        }
        i0++;
    }
}