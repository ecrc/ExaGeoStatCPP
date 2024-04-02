
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file UnivariateMaternDdsigmaSquare.cpp
 * @brief Implementation of the UnivariateMaternDdsigmaSquare kernel.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-04-14
**/

#include <kernels/concrete/UnivariateMaternDdsigmaSquare.hpp>

using namespace exageostat::kernels;
using namespace exageostat::dataunits;

template<typename T>
UnivariateMaternDdsigmaSquare<T>::UnivariateMaternDdsigmaSquare() {
    this->mP = 1;
    this->mParametersNumber = 3;
}


template<typename T>
Kernel<T> *UnivariateMaternDdsigmaSquare<T>::Create() {
    KernelsConfigurations::GetParametersNumberKernelMap()["UnivariateMaternDdsigmaSquare"] = 3;
    return new UnivariateMaternDdsigmaSquare();
}

namespace exageostat::kernels {
    template<typename T> bool UnivariateMaternDdsigmaSquare<T>::plugin_name = plugins::PluginRegistry<exageostat::kernels::Kernel<T>>::Add(
            "UnivariateMaternDdsigmaSquare", UnivariateMaternDdsigmaSquare::Create);
}

template<typename T>
void UnivariateMaternDdsigmaSquare<T>::GenerateCovarianceMatrix(T *apMatrixA, const int &aRowsNumber,
                                                                const int &aColumnsNumber, const int &aRowOffset,
                                                                const int &aColumnOffset, Locations<T> &aLocation1,
                                                                Locations<T> &aLocation2, Locations<T> &aLocation3,
                                                                T *aLocalTheta, const int &aDistanceMetric) {
    int i, j;
    int i0 = aRowOffset;
    int j0;
    T expr;
    T con;
    con = pow(2, (aLocalTheta[2] - 1)) * tgamma(aLocalTheta[2]);
    con = 1.0 / con;

    for (i = 0; i < aRowsNumber; i++) {
        j0 = aColumnOffset;
        for (j = 0; j < aColumnsNumber; j++) {
            expr = helpers::DistanceCalculationHelpers<T>::CalculateDistance(aLocation1, aLocation2, i0, j0,
                                                                             aDistanceMetric, 0) / aLocalTheta[1];
            if (expr == 0) {
                *(apMatrixA + i + j * aRowsNumber) = 1;
            } else {

                *(apMatrixA + i + j * aRowsNumber) = con * pow(expr, aLocalTheta[2]) * gsl_sf_bessel_Knu(aLocalTheta[2],
                                                                                                         expr); // derivative with respect to sigma square

            }
            j0++;

        }
        i0++;
    }
}
