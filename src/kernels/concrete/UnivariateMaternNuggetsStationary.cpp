
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file UnivariateMaternNuggetsStationary.cpp
 * @brief Implementation of the UnivariateMaternNuggetsStationary kernel.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-04-14
**/

#include <kernels/concrete/UnivariateMaternNuggetsStationary.hpp>


using namespace exageostat::kernels;
using namespace exageostat::dataunits;
using namespace exageostat::helpers;

template<typename T>
UnivariateMaternNuggetsStationary<T>::UnivariateMaternNuggetsStationary() {
    this->mP = 1;
    this->mParametersNumber = 4;
}

template<typename T>
Kernel<T> *UnivariateMaternNuggetsStationary<T>::Create() {
    KernelsConfigurations::GetParametersNumberKernelMap()["UnivariateMaternNuggetsStationary"] = 4;
    return new UnivariateMaternNuggetsStationary();
}

namespace exageostat::kernels {
    template<typename T> bool UnivariateMaternNuggetsStationary<T>::plugin_name = plugins::PluginRegistry<exageostat::kernels::Kernel<T>>::Add(
            "UnivariateMaternNuggetsStationary", UnivariateMaternNuggetsStationary<T>::Create);
}

template<typename T>
void UnivariateMaternNuggetsStationary<T>::GenerateCovarianceMatrix(T *apMatrixA, const int &aRowsNumber,
                                                                    const int &aColumnsNumber, const int &aRowOffset,
                                                                    const int &aColumnOffset, Locations<T> &aLocation1,
                                                                    Locations<T> &aLocation2, Locations<T> &aLocation3,
                                                                    T *aLocalTheta, const int &aDistanceMetric) {

    int i, j;
    int i0 = aRowOffset;
    int j0;
    T expr;
    T con;
    T sigma_square = aLocalTheta[0];

    con = pow(2, (aLocalTheta[2] - 1)) * tgamma(aLocalTheta[2]);
    con = 1.0 / con;
    con = sigma_square * con;
    int flag = aLocation1.GetLocationZ() == nullptr ? 0 : 1;

    if (aLocation1.GetLocationZ() == nullptr || aLocation2.GetLocationZ() == nullptr) {
        for (i = 0; i < aRowsNumber; i++) {
            j0 = aColumnOffset;
            for (j = 0; j < aColumnsNumber; j++) {
                expr = DistanceCalculationHelpers<T>::CalculateDistance(aLocation1, aLocation2, j0, i0, aDistanceMetric,
                                                                        flag) / aLocalTheta[1];
                if (expr == 0) {
                    apMatrixA[i + j * aRowsNumber] = sigma_square + aLocalTheta[3];
                } else {
                    // Matern Function
                    apMatrixA[i + j * aRowsNumber] =
                            con * pow(expr, aLocalTheta[2]) * gsl_sf_bessel_Knu(aLocalTheta[2], expr);
                }
                j0++;
            }
            i0++;
        }
    } else {
        for (i = 0; i < aRowsNumber; i++) {
            j0 = aColumnOffset;
            for (j = 0; j < aColumnsNumber; j++) {
                flag = 1;
                expr = DistanceCalculationHelpers<T>::CalculateDistance(aLocation1, aLocation2, j0, i0, aDistanceMetric,
                                                                        flag);
                if (expr == 0) {
                    apMatrixA[i + j * aRowsNumber] = sigma_square + aLocalTheta[3];
                } else {
                    // Matern Function
                    apMatrixA[i + j * aRowsNumber] =
                            con * pow(expr, aLocalTheta[2]) * gsl_sf_bessel_Knu(aLocalTheta[2], expr);
                }
                j0++;
            }
            i0++;
        }
    }
}