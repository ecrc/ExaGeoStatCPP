
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file DataGenerationFactory.hpp
 *
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-03-04
**/

#ifndef EXAGEOSTAT_CPP_DATAGENERATIONFACTORY_HPP
#define EXAGEOSTAT_CPP_DATAGENERATIONFACTORY_HPP

#include <data-generators/DataGenerator.hpp>
#include <configurations/data-generation/concrete/SyntheticDataConfigurations.h>

namespace exageostat {
    namespace generators {

        class DataGenerationFactory{
        public:
            DataGenerator *createGenerator(configurations::data_configurations::SyntheticDataConfigurations *aConfigurations);
        };

    }//namespace generators
}//namespace exageostat

#endif //EXAGEOSTAT_CPP_DATAGENERATIONFACTORY_HPP
