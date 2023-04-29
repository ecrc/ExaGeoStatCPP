// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file TestUnivariateMaternStationary.cpp
 *
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-04-29
**/
#include <libraries/catch/catch.hpp>
#include <linear-algebra-solvers/LinearAlgebraFactory.hpp>
#include <configurations/data-generation/concrete/SyntheticDataConfigurations.hpp>

using namespace exageostat::configurations::data_configurations;
using namespace exageostat::linearAlgebra;
using namespace exageostat::common;

void TEST_GENERATION(){
//    auto *syntheticDataConfigurations = new SyntheticDataConfigurations();
//#ifdef EXAGEOSTAT_USE_CHAMELEON
//    auto linearAlgebraSolver = LinearAlgebraFactory<double>::CreateLinearAlgebraSolver(EXACT_DENSE);
//#endif
//#ifdef EXAGEOSTAT_USE_HiCMA
//    auto linearAlgebraSolver = LinearAlgebraFactory<double>::CreateLinearAlgebraSolver(TILE_LOW_RANK);
//#endif
//    syntheticDataConfigurations->SetProblemSize(9);
//    syntheticDataConfigurations->SetDenseTileSize(5);
//    syntheticDataConfigurations->SetDimension(Dimension2D);
//    syntheticDataConfigurations->SetIsSynthetic(true);
//
//    linearAlgebraSolver->SetConfigurations(syntheticDataConfigurations);
//
//    linearAlgebraSolver->InitiateDescriptors();
//    linearAlgebraSolver->testKernelfornow();
}
TEST_CASE("Univariate Matern Stationary kernel test") {

    TEST_GENERATION();
}