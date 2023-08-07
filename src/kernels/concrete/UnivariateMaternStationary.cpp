
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file UnivariateMaternStationary.cpp
 *
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-04-14
**/

#include <kernels/concrete/UnivariateMaternStationary.hpp>

using namespace exageostat::kernels;
using namespace exageostat::dataunits;
using namespace std;

template<typename T>
UnivariateMaternStationary<T>::UnivariateMaternStationary() {
    this->mP = 1;
    this->mParametersNumber = 3;
}


template<typename T>
Kernel<T> *UnivariateMaternStationary<T>::Create() {
    return new UnivariateMaternStationary();
}

namespace exageostat::kernels {
    template<typename T> bool UnivariateMaternStationary<T>::plugin_name = plugins::PluginRegistry<exageostat::kernels::Kernel<T>>::Add(
            "UnivariateMaternStationary", UnivariateMaternStationary<T>::Create);
}

template<typename T>
void UnivariateMaternStationary<T>::GenerateCovarianceMatrix(T *apMatrixA, int &aRowsNumber, int &aColumnsNumber,
                                                          int &aRowOffset, int &aColumnOffset, Locations<T> *apLocation1,
                                                          Locations<T> *apLocation2, Locations<T> *apLocation3,
                                                          T *aLocalTheta, int &aDistanceMetric) {
    const T sigma_square = aLocalTheta[0];
    const T nu = aLocalTheta[2];
    const T inv_con = sigma_square * (1.0/(pow(2, (nu - 1)) * tgamma((nu))));

    int i0 = aRowOffset;
    int flag = 0;

    for (int i = 0; i < aRowsNumber; i++) {
        int j0 = aColumnOffset;
        for (int j = 0; j < aColumnsNumber; j++) {
            const T dist = this->CalculateDistance(apLocation1, apLocation2, i0, j0, aDistanceMetric, flag) / aLocalTheta[1];
            *(apMatrixA + i + j * aRowsNumber) = (dist == 0.0)
                                                 ? sigma_square
                                                 : inv_con * pow(dist, nu) * gsl_sf_bessel_Knu(nu, dist);
            j0++;
        }
        i0++;
    }
}
