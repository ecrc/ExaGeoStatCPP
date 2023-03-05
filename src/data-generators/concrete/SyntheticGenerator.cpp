
/**
 * @file SyntheticGenerator.cpp
 *
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-02-14
**/

#include <data-generators/concrete/SyntheticGenerator.hpp>

using namespace exageostat::generators::Synthetic;
using namespace exageostat::dataunits;
using namespace exageostat::configurations::data_configurations;

Locations SyntheticGenerator::InitializeLocations(Locations aLocations) {

    return aLocations;
}

void SyntheticGenerator::GenerateLocations(int aN, int aSeed) {

}

void SyntheticGenerator::Print() {
    std::cout << "HELLO YOU'RE USING SYNTHETIC DATA GENERATION" << std::endl;
}

SyntheticGenerator::SyntheticGenerator(SyntheticDataConfigurations *aConfigurations) {

    std::cout << "Hi: " << aConfigurations->GetProblemSize() << std::endl;
}
