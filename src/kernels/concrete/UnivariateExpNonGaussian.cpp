
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file UnivariateExpNonGaussian.cpp
 * @brief Implementation of the UnivariateExpNonGaussian kernel.
 * @version 1.0.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-04-14
**/

#include <kernels/concrete/UnivariateExpNonGaussian.hpp>

using namespace std;

using namespace exageostat::kernels;
using namespace exageostat::dataunits;

template<typename T>
UnivariateExpNonGaussian<T>::UnivariateExpNonGaussian() {
    this->mP = 1;
    this->mParametersNumber = 6;
}

template<typename T>
Kernel<T> *UnivariateExpNonGaussian<T>::Create() {
    return new UnivariateExpNonGaussian();
}

namespace exageostat::kernels {
    template<typename T> bool UnivariateExpNonGaussian<T>::plugin_name = plugins::PluginRegistry<exageostat::kernels::Kernel<T>>::Add(
            "UnivariateExpNonGaussian", UnivariateExpNonGaussian<T>::Create);
}

template<typename T>
void
UnivariateExpNonGaussian<T>::GenerateCovarianceMatrix(T *apMatrixA, const int &aRowsNumber, const int &aColumnsNumber,
                                                      const int &aRowOffset, const int &aColumnOffset,
                                                      dataunits::Locations<T> &aLocation1,
                                                      dataunits::Locations<T> &aLocation2,
                                                      dataunits::Locations<T> &aLocation3, T *aLocalTheta,
                                                      const int &aDistanceMetric) {

    int i, j;
    int i0 = aRowOffset;
    int j0;
    double expr;
    double sigma_square = 1;
    int flag = 0;

    for (i = 0; i < aRowsNumber; i++) {
        j0 = aColumnOffset;
        for (j = 0; j < aColumnsNumber; j++) {
            expr = this->CalculateDistance(aLocation1, aLocation2, i0, j0, aDistanceMetric, flag) / aLocalTheta[0];

            if (expr == 0)
                apMatrixA[i + j * aRowsNumber] = sigma_square /*+ 1e-4*/;
            else
                apMatrixA[i + j * aRowsNumber] = exp(-expr);

            j0++;
        }
        i0++;
    }
}