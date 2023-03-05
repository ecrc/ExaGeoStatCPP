
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file DataGenerationFactory.cpp
 *
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-03-04
**/
#include <data-generators/DataGenerationFactory.hpp>
#include <data-generators/concrete/SyntheticGenerator.hpp>

using namespace exageostat::generators;
using namespace exageostat::configurations::data_configurations;
using namespace exageostat::generators::Synthetic;

DataGenerator *DataGenerationFactory::createGenerator(SyntheticDataConfigurations *aConfigurations) {

    bool isSynthetic = aConfigurations->GetIsSynthetic();

    if (isSynthetic){
        return new SyntheticGenerator;
    }

    return nullptr;
}
