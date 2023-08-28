
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file UnivariateMaternDdsigmaSquare.cpp
 * @brief Implementation of the UnivariateMaternDdsigmaSquare kernel.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @author Mahmoud ElKarargy
 * @date 2023-04-14
**/

#include <kernels/concrete/UnivariateMaternDdsigmaSquare.hpp>

using namespace std;

using namespace exageostat::kernels;
using namespace exageostat::dataunits;

template<typename T>
UnivariateMaternDdsigmaSquare<T>::UnivariateMaternDdsigmaSquare() {
    this->mP = 1;
    this->mParametersNumber = 3;
}


template<typename T>
Kernel<T> *UnivariateMaternDdsigmaSquare<T>::Create() {
    return new UnivariateMaternDdsigmaSquare();
}

namespace exageostat::kernels {
    template<typename T> bool UnivariateMaternDdsigmaSquare<T>::plugin_name = plugins::PluginRegistry<exageostat::kernels::Kernel<T>>::Add(
            "UnivariateMaternDdsigmaSquare", UnivariateMaternDdsigmaSquare::Create);
}

template<typename T>
void UnivariateMaternDdsigmaSquare<T>::GenerateCovarianceMatrix(T *apMatrixA, const int &aRowsNumber,
                                                                const int &aColumnsNumber,
                                                                const int &aRowOffset, const int &aColumnOffset,
                                                                dataunits::Locations<T> &aLocation1,
                                                                dataunits::Locations<T> &aLocation2,
                                                                dataunits::Locations<T> &aLocation3, T *aLocalTheta,
                                                                const int &aDistanceMetric) {

    int i, j;
    //// TODO: Implementation is Empty in the old version!
    for (i = 0; i < aRowsNumber; i++) {
        for (j = 0; j < aColumnsNumber; j++) {
            apMatrixA[i + j * aRowsNumber] = 0.0;
        }
    }
}
