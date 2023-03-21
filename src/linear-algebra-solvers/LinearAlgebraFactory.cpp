
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file LinearAlgebraFactory.cpp
 *
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-03-20
**/

#include <linear-algebra-solvers/LinearAlgebraFactory.hpp>

using namespace exageostat::linearAlgebra;
using namespace exageostat::common;
using namespace exageostat::configurations;

template<typename T> std::unique_ptr<LinearAlgebraMethods<T>> LinearAlgebraFactory<T>::CreateLinearAlgebraSolver(Computation aComputation) {

    // Check the used Linear Algebra solver library, Whether it's HiCMA or Chameleon.

    // Chameleon Used
    if (aComputation == EXACT_DENSE) {
#ifdef EXAGEOSTAT_USE_CHAMELEON
        return std::make_unique<dense::ChameleonImplementation<T>>();
#else
        throw std::runtime_error("Dense matrix generation isn't supported without enabling Chameleon. Use -DEXAGEOSTAT_USE_CHAMELEON=ON");
        return nullptr;
#endif
    }

    // Hicma Used
    else if (aComputation == TILE_LOW_RANK) {
//        return std::make_unique<HiCMAAllocateDescriptors>();
    }
    //// TODO: which is Used?
    else if (aComputation == DIAGONAL_APPROX) {

    }

    // Return DataGenerator unique pointer of real type
    return nullptr;
}
