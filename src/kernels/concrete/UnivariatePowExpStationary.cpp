
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file UnivariatePowExpStationary.cpp
 * @brief Implementation of the UnivariatePowExpStationary kernel.
 * @version 1.1.0
 * @author Sameh Abdulah
 * @author Mahmoud ElKarargy
 * @date 2023-04-14
**/

#include <kernels/concrete/UnivariatePowExpStationary.hpp>
#include <helpers/DistanceCalculationHelpers.hpp>

using namespace exageostat::kernels;
using namespace exageostat::dataunits;
using namespace exageostat::helpers;

template<typename T>
UnivariatePowExpStationary<T>::UnivariatePowExpStationary() {
    this->mP = 1;
    this->mParametersNumber = 3;
}

template<typename T>
Kernel<T> *UnivariatePowExpStationary<T>::Create() {
    KernelsConfigurations::GetParametersNumberKernelMap()["UnivariatePowExpStationary"] = 3;
    return new UnivariatePowExpStationary();
}

namespace exageostat::kernels {
    template<typename T> bool UnivariatePowExpStationary<T>::plugin_name = plugins::PluginRegistry<exageostat::kernels::Kernel<T>>::Add(
            "UnivariatePowExpStationary", UnivariatePowExpStationary<T>::Create);
}

template<typename T>
void
UnivariatePowExpStationary<T>::GenerateCovarianceMatrix(T *apMatrixA, const int &aRowsNumber, const int &aColumnsNumber,
                                                        const int &aRowOffset, const int &aColumnOffset,
                                                        Locations<T> &aLocation1, Locations<T> &aLocation2,
                                                        Locations<T> &aLocation3, T *aLocalTheta,
                                                        const int &aDistanceMetric) {

    const T sigma_square = aLocalTheta[0];
    const T nu = aLocalTheta[2];
    int i0 = aRowOffset;
    int flag = 0;
    int j0;
    int i, j;
    T dist;

    for (i = 0; i < aRowsNumber; i++) {
        j0 = aColumnOffset;
        for (j = 0; j < aColumnsNumber; j++) {
            dist = DistanceCalculationHelpers<T>::CalculateDistance(aLocation1, aLocation2, i0, j0, aDistanceMetric,
                                                                    flag);
            dist = pow(dist, nu);
            *(apMatrixA + i + j * aRowsNumber) = (dist == 0.0) ? sigma_square : sigma_square *
                                                                                exp(-(dist / aLocalTheta[1]));
            j0++;
        }
        i0++;
    }
}
