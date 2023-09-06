
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file UnivariateMaternDdnuNu.cpp
 * @brief Implementation of the UnivariateMaternDdnuNu kernel.
 * @version 1.0.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-04-14
**/

#include <kernels/concrete/UnivariateMaternDdnuNu.hpp>

using namespace exageostat::kernels;
using namespace exageostat::dataunits;

template<typename T>
UnivariateMaternDdnuNu<T>::UnivariateMaternDdnuNu() {
    this->mP = 1;
    this->mParametersNumber = 3;
}

template<typename T>
Kernel<T> *UnivariateMaternDdnuNu<T>::Create() {
    KernelsConfigurations::GetParametersNumberKernelMap()["UnivariateMaternDdnuNu"] = 3;
    return new UnivariateMaternDdnuNu();
}

namespace exageostat::kernels {
    template<typename T> bool UnivariateMaternDdnuNu<T>::plugin_name = plugins::PluginRegistry<exageostat::kernels::Kernel<T>>::Add(
            "UnivariateMaternDdnuNu", UnivariateMaternDdnuNu<T>::Create);
}

template<typename T>
void
UnivariateMaternDdnuNu<T>::GenerateCovarianceMatrix(T *apMatrixA, const int &aRowsNumber, const int &aColumnsNumber,
                                                    const int &aRowOffset, const int &aColumnOffset,
                                                    Locations<T> &aLocation1, Locations<T> &aLocation2,
                                                    Locations<T> &aLocation3, T *aLocalTheta,
                                                    const int &aDistanceMetric) {

    int i, j;
    int i0 = aRowOffset;
    int j0;
    double expr;
    double con;
    double nu_expr;
    double nu_expr_dprime;
    double sigma_square = aLocalTheta[0];
    con = pow(2, (aLocalTheta[2] - 1)) * tgamma(aLocalTheta[2]);
    con = 1.0 / con;
    int flag = 0;

    for (i = 0; i < aRowsNumber; i++) {
        j0 = aColumnOffset;
        for (j = 0; j < aColumnsNumber; j++) {
            expr = this->CalculateDistance(aLocation1, aLocation2, i0, j0, aDistanceMetric, flag) / aLocalTheta[1];
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
                                                                      this->CalculateDerivativeBesselNu(aLocalTheta[2],
                                                                                                        expr)));
                nu_expr_dprime = (1 - aLocalTheta[2]) * 1 / pow(2, aLocalTheta[2]) * 1 / tgamma(aLocalTheta[2]) *
                                 pow(expr, aLocalTheta[2]) * this->CalculateDerivativeBesselNu(aLocalTheta[2], expr) +
                                 pow(2, 1 - aLocalTheta[2]) *
                                 (-1 / tgamma(aLocalTheta[2]) * gsl_sf_psi(aLocalTheta[2]) * pow(expr, aLocalTheta[2]) *
                                  this->CalculateDerivativeBesselNu(aLocalTheta[2], expr) + 1 / tgamma(aLocalTheta[2]) *
                                                                                            (pow(expr, aLocalTheta[2]) *
                                                                                             log(expr) *
                                                                                             this->CalculateDerivativeBesselNu(
                                                                                                     aLocalTheta[2],
                                                                                                     expr) +
                                                                                             pow(expr, aLocalTheta[2]) *
                                                                                             this->CalculateSecondDerivativeBesselNu(
                                                                                                     aLocalTheta[2],
                                                                                                     expr)));
                apMatrixA[i + j * aRowsNumber] =
                        (-0.5 * con * pow(expr, aLocalTheta[2]) * gsl_sf_bessel_Knu(aLocalTheta[2], expr) +
                         (1 - aLocalTheta[2]) / 2 * nu_expr -
                         (1 / aLocalTheta[2] + 0.5 * pow(aLocalTheta[2], 2)) * con * pow(expr, aLocalTheta[2]) *
                         gsl_sf_bessel_Knu(aLocalTheta[2], expr) - gsl_sf_psi(aLocalTheta[2]) * nu_expr +
                         log(expr) * nu_expr + nu_expr_dprime) * sigma_square;
            }
            j0++;
        }
        i0++;
    }
}