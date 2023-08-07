
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
using namespace exageostat::generators::synthetic;
using namespace exageostat::dataunits;
using namespace exageostat::configurations;

template<typename T>
std::unique_ptr<DataGenerator<T>> DataGenerator<T>::CreateGenerator(Configurations *apConfigurations) {

    // Check the used Data generation method, whether it's synthetic or real.
    bool is_synthetic = apConfigurations->GetIsSynthetic();

    // Return DataGenerator unique pointer of Synthetic type
    if (is_synthetic) {
        return std::unique_ptr<DataGenerator<T>>(SyntheticGenerator<T>::GetInstance(apConfigurations));
    } else {
        // Open saved files
        throw std::domain_error("Unsupported for now, Please add --synthetic_data");
    }
}

template<typename T>
Locations<T> *DataGenerator<T>::GetLocations() {
    if (this->mpLocations == nullptr) {
        throw std::runtime_error("Locations is null");
    }
    return this->mpLocations;
}

template<typename T>
exageostat::kernels::Kernel<T> *DataGenerator<T>::GetKernel() {
    if (this->mpKernel == nullptr) {
        throw std::runtime_error("Kernel is null");
    }
    return this->mpKernel;
}

template<typename T>
DataGenerator<T>::~DataGenerator(){

    // Check the used Data generation method, whether it's synthetic or real.
    bool is_synthetic = this->mpConfigurations->GetIsSynthetic();

    // Return DataGenerator unique pointer of Synthetic type
    if (is_synthetic) {
        SyntheticGenerator<T>::GetInstance(this->mpConfigurations)->ReleaseInstance();
    } else {
        // Open saved files
        throw std::domain_error("Unsupported for now, Please add --synthetic_data");
    }
}