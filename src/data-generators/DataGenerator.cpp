
/**
 * @file DataGenerator.cpp
 * @brief Implementation of DataGenerator class
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-02-14
**/

#include <data-generators/DataGenerator.hpp>
#include <data-generators/concrete/SyntheticGenerator.hpp>

using namespace exageostat::generators;
using namespace exageostat::configurations::data_configurations;
using namespace exageostat::generators::Synthetic;
using namespace exageostat::dataunits;

std::unique_ptr<DataGenerator> DataGenerator::CreateGenerator(SyntheticDataConfigurations *apConfigurations) {
    // Check the used Data generation method, whether it's synthetic or real.
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

void DataGenerator::InitLocationsClass() {
    this->mpLocations = new Locations();
}

Locations* DataGenerator::GetLocations() {
    return this->mpLocations;
}