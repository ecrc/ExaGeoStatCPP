
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file UnivariateMaternNonGaussian.cpp
 * @brief Implementation of the UnivariateMaternNonGaussian kernel.
 * @version 1.0.1
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-04-14
**/

#include <kernels/concrete/UnivariateMaternNonGaussian.hpp>
#include <helpers/DistanceCalculationHelpers.hpp>

using namespace exageostat::kernels;
using namespace exageostat::dataunits;
using namespace exageostat::helpers;

template<typename T>
UnivariateMaternNonGaussian<T>::UnivariateMaternNonGaussian() {
    this->mP = 1;
    this->mParametersNumber = 6;
}

template<typename T>
Kernel<T> *UnivariateMaternNonGaussian<T>::Create() {
    KernelsConfigurations::GetParametersNumberKernelMap()["UnivariateMaternNonGaussian"] = 6;
    return new UnivariateMaternNonGaussian();
}

namespace exageostat::kernels {
    template<typename T> bool UnivariateMaternNonGaussian<T>::plugin_name = plugins::PluginRegistry<exageostat::kernels::Kernel<T>>::Add(
            "UnivariateMaternNonGaussian", UnivariateMaternNonGaussian<T>::Create);
}

template<typename T>
void UnivariateMaternNonGaussian<T>::GenerateCovarianceMatrix(T *apMatrixA, const int &aRowsNumber,
                                                              const int &aColumnsNumber, const int &aRowOffset,
                                                              const int &aColumnOffset, Locations<T> &aLocation1,
                                                              Locations<T> &aLocation2, Locations<T> &aLocation3,
                                                              T *aLocalTheta, const int &aDistanceMetric) {

    int i, j;
    int i0 = aRowOffset;
    int j0;
    T expr;
    T con;
    T sigma_square = 1;

    con = pow(2, (aLocalTheta[1] - 1)) * tgamma(aLocalTheta[1]);
    con = 1.0 / con;
    con = sigma_square * con;
    int flag = aLocation1.GetLocationZ() == nullptr ? 0 : 1;

    for (i = 0; i < aRowsNumber; i++) {
        j0 = aColumnOffset;
        for (j = 0; j < aColumnsNumber; j++) {
            expr = 4 * sqrt(2 * aLocalTheta[1]) *
                   (DistanceCalculationHelpers<T>::CalculateDistance(aLocation1, aLocation2, i0, j0, aDistanceMetric,
                                                                     flag) / aLocalTheta[0]);
            if (expr == 0) {
                apMatrixA[i + j * aRowsNumber] = sigma_square;
            } else {
                // Matern Function
                apMatrixA[i + j * aRowsNumber] =
                        con * pow(expr, aLocalTheta[1]) * gsl_sf_bessel_Knu(aLocalTheta[1], expr);
            }
            j0++;
        }
        i0++;
    }
}