
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file UnivariateMaternDdsigmaSquareNu.cpp
 * @brief Implementation of the UnivariateMaternDdsigmaSquareNu kernel.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @author Mahmoud ElKarargy
 * @date 2023-04-14
**/

#include <kernels/concrete/UnivariateMaternDdsigmaSquareNu.hpp>

using namespace std;

using namespace exageostat::kernels;
using namespace exageostat::dataunits;

template<typename T>
UnivariateMaternDdsigmaSquareNu<T>::UnivariateMaternDdsigmaSquareNu() {
    this->mP = 1;
    this->mParametersNumber = 3;
}

template<typename T>
Kernel<T> *UnivariateMaternDdsigmaSquareNu<T>::Create() {
    return new UnivariateMaternDdsigmaSquareNu();
}

namespace exageostat::kernels {
    template<typename T> bool UnivariateMaternDdsigmaSquareNu<T>::plugin_name = plugins::PluginRegistry<exageostat::kernels::Kernel<T>>::Add(
            "UnivariateMaternDdsigmaSquareNu", UnivariateMaternDdsigmaSquareNu<T>::Create);
}

template<typename T>
void UnivariateMaternDdsigmaSquareNu<T>::GenerateCovarianceMatrix(T *apMatrixA, int &aRowsNumber, int &aColumnsNumber,
                                                                  int &aRowOffset, int &aColumnOffset,
                                                                  Locations<T> *apLocation1,
                                                                  Locations<T> *apLocation2, Locations<T> *apLocation3,
                                                                  T *aLocalTheta, int &aDistanceMetric) {

    int i, j;
    int i0 = aRowOffset;
    int j0;
    double expr;
    double nu_expr;
    int flag = 0;

    for (i = 0; i < aRowsNumber; i++) {
        j0 = aColumnOffset;
        for (j = 0; j < aColumnsNumber; j++) {
            expr = this->CalculateDistance(*apLocation1, *apLocation2, i0, j0, aDistanceMetric, flag) / aLocalTheta[1];
            if (expr == 0) {
                apMatrixA[i + j * aRowsNumber] = 0.0;
            } else {
                nu_expr = (1 - aLocalTheta[2]) * 1 / pow(2, aLocalTheta[2]) * 1 / tgamma(aLocalTheta[2])
                          * pow(expr, aLocalTheta[2]) * gsl_sf_bessel_Knu(aLocalTheta[2], expr) +
                          pow(2, 1 - aLocalTheta[2])
                          * (-1 / tgamma(aLocalTheta[2]) * gsl_sf_psi(aLocalTheta[2]) * pow(expr, aLocalTheta[2])
                             * gsl_sf_bessel_Knu(aLocalTheta[2], expr) + 1 / tgamma(aLocalTheta[2])
                                                                         * (pow(expr, aLocalTheta[2]) * log(expr)
                                                                            * gsl_sf_bessel_Knu(aLocalTheta[2], expr) +
                                                                            pow(expr, aLocalTheta[2])
                                                                            * (gsl_sf_bessel_Knu(
                                                                                    aLocalTheta[2] + 0.000000001,
                                                                                    expr) -
                                                                               gsl_sf_bessel_Knu(aLocalTheta[2],
                                                                                                 expr)) / 0.000000001));
                apMatrixA[i + j * aRowsNumber] = nu_expr; // derivative with respect to sigma square
            }
            j0++;
        }
        i0++;
    }
}