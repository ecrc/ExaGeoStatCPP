
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ChameleonImplementationDense.cpp
 * @brief Dense Tile implementation of linear algebra methods.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-03-20
**/

#include <mkl_service.h>

#include <linear-algebra-solvers/concrete/chameleon/dense/ChameleonDense.hpp>
#include <linear-algebra-solvers/concrete/ChameleonHeaders.hpp>

using namespace std;

using namespace exageostat::linearAlgebra::dense;

template<typename T>
void
ChameleonDense<T>::ExaGeoStatPotrfTile(const common::UpperLower &aUpperLower, void *apA, int aBand,
                                       void *apCD, void *apCrk, const int &aMaxRank, const int &aAcc) {
    int status = CHAMELEON_dpotrf_Tile((cham_uplo_t) aUpperLower, (CHAM_desc_t *) apA);
    if (status != CHAMELEON_SUCCESS) {
        throw std::runtime_error("CHAMELEON_dpotrf_Tile Failed, Matrix is not positive definite");
    }
    // Due to a leak in dense mode in Chameleon, We had to free the buffer manually.
    mkl_free_buffers();
}