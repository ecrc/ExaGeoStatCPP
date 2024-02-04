
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file UnivariateMaternStationary.cpp
 * @brief Implementation of the UnivariateMaternStationary kernel.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-04-14
**/

#include <kernels/concrete/UnivariateMaternStationary.hpp>
#include <helpers/DistanceCalculationHelpers.hpp>

using namespace exageostat::kernels;
using namespace exageostat::dataunits;
using namespace exageostat::helpers;

template<typename T>
UnivariateMaternStationary<T>::UnivariateMaternStationary() {
    this->mP = 1;
    this->mParametersNumber = 3;
}

template<typename T>
Kernel<T> *UnivariateMaternStationary<T>::Create() {
    KernelsConfigurations::GetParametersNumberKernelMap()["UnivariateMaternStationary"] = 3;
    return new UnivariateMaternStationary();
}

namespace exageostat::kernels {
    template<typename T> bool UnivariateMaternStationary<T>::plugin_name = plugins::PluginRegistry<exageostat::kernels::Kernel<T>>::Add(
            "UnivariateMaternStationary", UnivariateMaternStationary<T>::Create);
}

template<typename T>
void
UnivariateMaternStationary<T>::GenerateCovarianceMatrix(T *apMatrixA, const int &aRowsNumber, const int &aColumnsNumber,
                                                        const int &aRowOffset, const int &aColumnOffset,
                                                        Locations<T> &aLocation1, Locations<T> &aLocation2,
                                                        Locations<T> &aLocation3, T *aLocalTheta,
                                                        const int &aDistanceMetric) {

    const T sigma_square = aLocalTheta[0];
    const T nu = aLocalTheta[2];
    const T inv_con = sigma_square * (1.0 / (pow(2, (nu - 1)) * tgamma((nu))));
    int i0 = aRowOffset;
    int flag = aLocation1.GetLocationZ() == nullptr ? 0 : 1;
    int j0;
    int i, j;
    T dist;

    for (i = 0; i < aRowsNumber; i++) {
        j0 = aColumnOffset;
        for (j = 0; j < aColumnsNumber; j++) {
            dist = DistanceCalculationHelpers<T>::CalculateDistance(aLocation1, aLocation2, i0, j0, aDistanceMetric,
                                                                    flag) / aLocalTheta[1];
            *(apMatrixA + i + j * aRowsNumber) = (dist == 0.0) ? sigma_square : inv_con * pow(dist, nu) *
                                                                                gsl_sf_bessel_Knu(nu, dist);
            j0++;
        }
        i0++;
    }
}