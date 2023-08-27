
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file UnivariateMaternDsigmaSquare.cpp
 * @brief Implementation of the UnivariateMaternDsigmaSquare kernel.
  * @version 1.0.0
 * @author Sameh Abdulah* @author Suhas Shankar
 * @author Mary Lai Salvana
 * @author Mahmoud ElKarargy
 * @date 2023-04-14
**/

#include <kernels/concrete/UnivariateMaternDsigmaSquare.hpp>

using namespace std;

using namespace exageostat::kernels;
using namespace exageostat::dataunits;

template<typename T>
UnivariateMaternDsigmaSquare<T>::UnivariateMaternDsigmaSquare() {
    this->mP = 1;
    this->mParametersNumber = 3;
}

template<typename T>
Kernel<T> *UnivariateMaternDsigmaSquare<T>::Create() {
    return new UnivariateMaternDsigmaSquare();
}

namespace exageostat::kernels {
    template<typename T> bool UnivariateMaternDsigmaSquare<T>::plugin_name = plugins::PluginRegistry<exageostat::kernels::Kernel<T>>::Add(
            "UnivariateMaternDsigmaSquare", UnivariateMaternDsigmaSquare<T>::Create);
}

template<typename T>
void UnivariateMaternDsigmaSquare<T>::GenerateCovarianceMatrix(T *apMatrixA, const int &aRowsNumber,
                                                               const int &aColumnsNumber,
                                                               const int &aRowOffset, const int &aColumnOffset,
                                                               dataunits::Locations<T> &aLocation1,
                                                               dataunits::Locations<T> &aLocation2,
                                                               dataunits::Locations<T> &aLocation3, T *aLocalTheta,
                                                               const int &aDistanceMetric) {
    int i, j;
    int i0 = aRowOffset;
    int j0;
    double expr;
    double con;

    con = pow(2, (aLocalTheta[2] - 1)) * tgamma(aLocalTheta[2]);
    con = 1.0 / con;
    int flag = 0;
    for (i = 0; i < aRowsNumber; i++) {
        j0 = aColumnOffset;
        for (j = 0; j < aColumnsNumber; j++) {
            expr = CalculateDistance(aLocation1, aLocation2, i0, j0, aDistanceMetric, flag) / aLocalTheta[1];
            if (expr == 0) {
                apMatrixA[i + j * aRowsNumber] = 1;
            } else {
                apMatrixA[i + j * aRowsNumber] = con * pow(expr, aLocalTheta[2]) *
                                                 gsl_sf_bessel_Knu(aLocalTheta[2],
                                                                   expr); // derivative with respect to sigma square
            }
            j0++;
        }
        i0++;
    }
}