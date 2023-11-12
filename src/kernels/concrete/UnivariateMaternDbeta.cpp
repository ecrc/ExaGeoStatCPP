
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file UnivariateMaternDbeta.cpp
 * @brief Implementation of the UnivariateMaternDbeta kernel.
 * @version 1.0.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-04-14
**/

#include <kernels/concrete/UnivariateMaternDbeta.hpp>
#include <helpers/DistanceCalculationHelpers.hpp>

using namespace exageostat::kernels;
using namespace exageostat::dataunits;
using namespace exageostat::helpers;

template<typename T>
UnivariateMaternDbeta<T>::UnivariateMaternDbeta() {
    this->mP = 1;
    this->mParametersNumber = 3;
}

template<typename T>
Kernel<T> *UnivariateMaternDbeta<T>::Create() {
    KernelsConfigurations::GetParametersNumberKernelMap()["UnivariateMaternDbeta"] = 3;
    return new UnivariateMaternDbeta();
}

namespace exageostat::kernels {
    template<typename T> bool UnivariateMaternDbeta<T>::plugin_name = plugins::PluginRegistry<exageostat::kernels::Kernel<T>>::Add(
            "UnivariateMaternDbeta", UnivariateMaternDbeta<T>::Create);
}

template<typename T>
void UnivariateMaternDbeta<T>::GenerateCovarianceMatrix(T *apMatrixA, const int &aRowsNumber, const int &aColumnsNumber,
                                                        const int &aRowOffset, const int &aColumnOffset,
                                                        Locations<T> &aLocation1, Locations<T> &aLocation2,
                                                        Locations<T> &aLocation3, T *aLocalTheta,
                                                        const int &aDistanceMetric) {

    int i, j;
    int i0 = aRowOffset;
    int j0;
    double expr;
    double beta_expr;
    double con;
    double sigma_square = aLocalTheta[0];
    con = pow(2, (aLocalTheta[2] - 1)) * tgamma(aLocalTheta[2]);
    con = 1.0 / con;
    int flag = 0;

    for (i = 0; i < aRowsNumber; i++) {
        j0 = aColumnOffset;
        for (j = 0; j < aColumnsNumber; j++) {
            expr = DistanceCalculationHelpers<T>::CalculateDistance(aLocation1, aLocation2, i0, j0, aDistanceMetric,
                                                                    flag) / aLocalTheta[1];
            if (expr == 0) {
                apMatrixA[i + j * aRowsNumber] = 0.0;
            } else {
                //derivative with respect to beta
                beta_expr = -aLocalTheta[2] / aLocalTheta[1] * pow(expr, aLocalTheta[2]) *
                            gsl_sf_bessel_Knu(aLocalTheta[2], expr) - pow(expr, aLocalTheta[2]) *
                                                                      (aLocalTheta[2] / expr *
                                                                       gsl_sf_bessel_Knu(aLocalTheta[2], expr) -
                                                                       gsl_sf_bessel_Knu(aLocalTheta[2] + 1, expr)) *
                                                                      expr / aLocalTheta[1];
                apMatrixA[i + j * aRowsNumber] = sigma_square * con * beta_expr;
            }
            j0++;
        }
        i0++;
    }
}