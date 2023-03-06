
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

std::unique_ptr<DataGenerator> DataGenerator::CreateGenerator(SyntheticDataConfigurations *apConfigurations) {

    // Check the used Data generation method, Whether it's synthetic or real.
    bool isSynthetic = apConfigurations->GetIsSynthetic();

    // Return DataGenerator unique pointer of Synthetic type
    if (isSynthetic) {
        return std::make_unique<SyntheticGenerator>(apConfigurations);
    }

    // Return DataGenerator unique pointer of real type
    return nullptr;
}

void DataGenerator::SetConfigurations(SyntheticDataConfigurations *apConfigurations) {
    this->mpConfigurations = apConfigurations;
}
