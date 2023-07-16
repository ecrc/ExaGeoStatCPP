
/*
 * Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
 * All rights reserved.
 * ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).
 */

/**
 * @file DataConfigurations.cpp
 * @brief This file contains the implementation of the DataConfigurations class.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-02-03
**/

#include <algorithm>
#include <vector>

#include <configurations/data-generation/DataConfigurations.hpp>

using namespace std;

using namespace exageostat::configurations::data_configurations;
using namespace exageostat::common;

void DataConfigurations::SetIsSynthetic(bool aIsSynthetic) {
    this->mIsSynthetic = aIsSynthetic;
}

bool DataConfigurations::GetIsSynthetic() const {
    return this->mIsSynthetic;
}

void DataConfigurations::SetLowerBounds(std::vector<double> &apTheta) {
    this->mLowerBounds = apTheta;
}

std::vector<double> &DataConfigurations::GetLowerBounds() {
    return this->mLowerBounds;
}

void DataConfigurations::SetUpperBounds(std::vector<double> &apTheta) {
    this->mUpperBounds = apTheta;
}

std::vector<double> &DataConfigurations::GetUpperBounds() {
    return this->mUpperBounds;
}

void DataConfigurations::SetStartingTheta(std::vector<double> &apTheta) {
    this->mStartingTheta = apTheta;
}

std::vector<double> &DataConfigurations::GetStartingTheta() {
    return this->mStartingTheta;
}

void DataConfigurations::SetTargetTheta(std::vector<double> &apTheta) {
    this->mTargetTheta = apTheta;
}

std::vector<double> &DataConfigurations::GetTargetTheta() {
    return this->mTargetTheta;
}

