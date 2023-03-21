
/*
 * Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
 * Copyright (C) 2023 by Brightskies inc,
 * All rights reserved.
 * ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).
 */

/**
 * @file DataConfigurations.cpp
 *
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-02-03
**/

#include <configurations/data-generation/DataConfigurations.hpp>

using namespace exageostat::configurations::data_configurations;
using namespace std;
using namespace exageostat::common;

string DataConfigurations::GetKernel() {
    return this->mKernel;
}

void DataConfigurations::SetKernel(std::string aKernel) {
    this->mKernel = aKernel;
}

void DataConfigurations::SetIsSynthetic(bool aIsSynthetic) {
    this->mIsSynthetic = aIsSynthetic;
}

bool DataConfigurations::GetIsSynthetic() {
    return this->mIsSynthetic;
}

void DataConfigurations::CheckKernelValue(std::string aKernel) {

    // finding position of input kernel
    auto position = availableKernels.find(aKernel);

    // If the element is not found,  then the iterator points to the position just after the last element in the set.
    if (position == availableKernels.end()) {
        throw range_error("Invalid value for Kernel. Please check manual.");
    }
}
