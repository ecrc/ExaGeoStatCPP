
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file StarPuHelpersFactory.cpp
 * @brief Factory for StarPu helpers.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @date 2024-02-25
**/

#include <runtime/starpu/helpers/StarPuHelpersFactory.hpp>
#include <runtime/starpu/helpers/concrete/ChameleonStarPuHelpers.hpp>

#ifdef USE_HICMA

#include <runtime/starpu/helpers/concrete/HicmaStarPuHelpers.hpp>

#endif

using namespace std;
using namespace exageostat::common;
using namespace exageostat::runtime;

unique_ptr<StarPuHelpers> StarPuHelpersFactory::CreateStarPuHelper(const Computation &aComputation) {
    if (aComputation == EXACT_DENSE || aComputation == DIAGONAL_APPROX) {
        return make_unique<ChameleonStarPuHelpers>();
    } else if (aComputation == TILE_LOW_RANK) {
#ifdef USE_HICMA
        return make_unique<HicmaStarPuHelpers>();
#else
        throw runtime_error("Tile low rank generation isn't supported without enabling HiCMA. Use -DUSE_HICMA=ON");
#endif
    }
    throw runtime_error("You need to enable whether HiCMA or Chameleon");
}
