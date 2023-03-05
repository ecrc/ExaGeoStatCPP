
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file SyntheticGenerator.cpp
 *
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-03-04
**/

#include <iostream>
#include <data-generators/DataGenerator.hpp>

using namespace exageostat::generators;
using namespace std;
using namespace exageostat::configurations::data_configurations;

int main(int argc, char **argv) {
    DataGenerator *dataGenerationFactory;
    DataGenerator *syntheticGenerator;

    SyntheticDataConfigurations *syntheticDataConfigurations = new SyntheticDataConfigurations(argc, argv);

    syntheticGenerator = dataGenerationFactory->CreateGenerator(syntheticDataConfigurations);
    syntheticGenerator->Print();

    return 0;
}