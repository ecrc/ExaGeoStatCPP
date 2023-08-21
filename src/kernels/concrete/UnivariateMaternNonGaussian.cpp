
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file UnivariateMaternNonGaussian.cpp
 * @brief Implementation of the UnivariateMaternNonGaussian kernel.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-04-14
**/

#include <kernels/concrete/UnivariateMaternNonGaussian.hpp>

using namespace std;

using namespace exageostat::kernels;
using namespace exageostat::dataunits;

template<typename T>
UnivariateMaternNonGaussian<T>::UnivariateMaternNonGaussian() {
    this->mP = 1;
    this->mParametersNumber = 6;
}

template<typename T>
Kernel<T> *UnivariateMaternNonGaussian<T>::Create() {
    return new UnivariateMaternNonGaussian();
}

namespace exageostat::kernels {
    template<typename T> bool UnivariateMaternNonGaussian<T>::plugin_name = plugins::PluginRegistry<exageostat::kernels::Kernel<T>>::Add(
            "UnivariateMaternNonGaussian", UnivariateMaternNonGaussian<T>::Create);
}

template<typename T>
void UnivariateMaternNonGaussian<T>::GenerateCovarianceMatrix(T *apMatrixA, int &aRowsNumber, int &aColumnsNumber,
                                                              int &aRowOffset, int &aColumnOffset,
                                                              Locations<T> *apLocation1,
                                                              Locations<T> *apLocation2, Locations<T> *apLocation3,
                                                              T *aLocalTheta, int &aDistanceMetric) {
    //localtheta[0] <- \phi
    //localtheta[1] <- \nu
    int i, j;
    int i0 = aRowOffset;
    int j0;
    double expr;
    double con;
    double sigma_square = 1;

    con = pow(2, (aLocalTheta[1] - 1)) * tgamma(aLocalTheta[1]);
    con = 1.0 / con;
    con = sigma_square * con;
    int flag = 0;

    for (i = 0; i < aRowsNumber; i++) {
        j0 = aColumnOffset;
        for (j = 0; j < aColumnsNumber; j++) {
            expr = 4 * sqrt(2 * aLocalTheta[1]) *
                   (this->CalculateDistance(*apLocation1, *apLocation2, i0, j0, aDistanceMetric, flag) /
                    aLocalTheta[0]);
            if (expr == 0)
                apMatrixA[i + j * aRowsNumber] = sigma_square /*+ 1e-4*/;
            else
                apMatrixA[i + j * aRowsNumber] = con * pow(expr, aLocalTheta[1])
                                                 * gsl_sf_bessel_Knu(aLocalTheta[1], expr); // Matern Function
            j0++;
        }
        i0++;
    }
}