
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

template<typename T>
std::unique_ptr<DataGenerator<T>> DataGenerator<T>::CreateGenerator(SyntheticDataConfigurations *apConfigurations) {
    // Check the used Data generation method, whether it's synthetic or real.
    bool is_synthetic = apConfigurations->GetIsSynthetic();

    // Return DataGenerator unique pointer of Synthetic type
    if (is_synthetic) {
        return std::make_unique<SyntheticGenerator<T>>(apConfigurations);
    } else {
        // Open saved files
        throw std::domain_error("Unsupported for now, Please add --SyntheticData");
    }
}

template<typename T>
Locations *DataGenerator<T>::GetLocations() {
    return this->mpLocations;
}

template<typename T>
exageostat::kernels::Kernel *DataGenerator<T>::GetKernel() {
    return this->mpKernel;
}
