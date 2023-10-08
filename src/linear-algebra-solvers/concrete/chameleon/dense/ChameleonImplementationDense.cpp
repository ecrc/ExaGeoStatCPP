
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ChameleonImplementationDense.cpp
 * @brief Dense Tile implementation of linear algebra methods.
 * @version 1.0.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-03-20
**/

#include <linear-algebra-solvers/concrete/chameleon/dense/ChameleonImplementationDense.hpp>

using namespace std;

using namespace exageostat::linearAlgebra::dense;
using namespace exageostat::common;
using namespace exageostat::dataunits;
using namespace exageostat::helpers;
using namespace exageostat::configurations;
using namespace exageostat::hardware;

template<typename T>
void
ChameleonImplementationDense<T>::ExaGeoStatPotrfTile(const common::UpperLower &aUpperLower, void *apA, int aDiagThick,
                                                     void *apCD, void *apCrk,
                                                     const int &aMaxRank, const int &aAcc) {
    int status = CHAMELEON_dpotrf_Tile((cham_uplo_t) aUpperLower, (CHAM_desc_t *) apA);
    if (status != CHAMELEON_SUCCESS) {
        throw std::runtime_error("CHAMELEON_dpotrf_Tile Failed, Matrix is not positive definite");
    }
}