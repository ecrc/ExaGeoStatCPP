
/**
 * @file DataGenerator.cpp
 *
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-02-14
**/

#include <data-generators/DataGenerator.hpp>
#include <data-generators/concrete/SyntheticGenerator.hpp>

using namespace exageostat::generators;
using namespace exageostat::configurations::data_configurations;
using namespace exageostat::generators::Synthetic;

DataGenerator *DataGenerator::createGenerator(SyntheticDataConfigurations *aConfigurations) {

    bool isSynthetic = aConfigurations->GetIsSynthetic();

    if (isSynthetic){
        return new SyntheticGenerator(aConfigurations);
    }

    return nullptr;
}
