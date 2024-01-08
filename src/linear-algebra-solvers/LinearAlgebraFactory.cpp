
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file LinearAlgebraFactory.cpp
 * @brief Implementation of the LinearAlgebraFactory class for creating linear algebra solvers for different computations using HiCMA or Chameleon libraries.
 * The factory creates a unique pointer to a concrete implementation of the LinearAlgebraMethods class based on the computation specified.
 * If the required library is not enabled, it throws a runtime_error exception.
 * @version 1.0.0
 * @author Mahmoud ElKarargy
 * @date 2023-03-20
**/

#include <linear-algebra-solvers/LinearAlgebraFactory.hpp>


#include <linear-algebra-solvers/concrete/chameleon/dense/ChameleonImplementationDense.hpp>
#include <linear-algebra-solvers/concrete/chameleon/diagonal-super-tile/ChameleonImplementationDST.hpp>

#ifdef USE_HICMA

#include <linear-algebra-solvers/concrete/hicma/tile-low-rank/HicmaImplementation.hpp>

#endif

using namespace exageostat::linearAlgebra;
using namespace exageostat::common;
using namespace exageostat::configurations;

template<typename T>
std::unique_ptr<LinearAlgebraMethods<T>> LinearAlgebraFactory<T>::CreateLinearAlgebraSolver(Computation aComputation) {

    // Check the used Linear Algebra solver library, whether it's HiCMA or Chameleon.
    if (aComputation == EXACT_DENSE) {

        return std::make_unique<dense::ChameleonImplementationDense<T>>();
    }

        // HiCMA Used
    else if (aComputation == TILE_LOW_RANK) {
#ifdef USE_HICMA
        return std::make_unique<tileLowRank::HicmaImplementation<T>>();
#else
        throw std::runtime_error(
                "Tile low rank generation isn't supported without enabling HiCMA. Use -DUSE_HICMA=ON");
#endif
    } else if (aComputation == DIAGONAL_APPROX) {
        return std::make_unique<diagonalSuperTile::ChameleonImplementationDST<T>>();

    }
    // Return nullptr if no computation is selected
    throw std::runtime_error("You need to enable whether HiCMA or Chameleon");
}
