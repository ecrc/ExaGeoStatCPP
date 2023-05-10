
/*
 * Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
 * Copyright (C) 2023 by Brightskies inc,
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

#include <configurations/data-generation/DataConfigurations.hpp>
#include <utility>
#include <algorithm>

using namespace exageostat::configurations::data_configurations;
using namespace std;
using namespace exageostat::common;

string DataConfigurations::GetKernel() {
    return this->mKernel;
}

void DataConfigurations::SetKernel(const std::string &aKernel) {
    this->mKernel = aKernel;
}

void DataConfigurations::SetIsSynthetic(bool aIsSynthetic) {
    this->mIsSynthetic = aIsSynthetic;
}

bool DataConfigurations::GetIsSynthetic() const {
    return this->mIsSynthetic;
}

void DataConfigurations::SetParametersNumber(int aParameterNumbers) {
    this->mParametersNumber = aParameterNumbers;
}

int DataConfigurations::GetParamtersNumber() {
    return this->mParametersNumber;
}

void DataConfigurations::SetLowerBounds(int aSize) {
    this->mLowerBounds =  (double *) std::malloc(aSize * sizeof(double ));
}

double *DataConfigurations::GetLowerBounds() {
    return this->mLowerBounds;
}

void DataConfigurations::SetUpperBounds(int aSize) {
    this->mUpperBounds =  (double *) std::malloc(aSize * sizeof(double ));
}

double *DataConfigurations::GetUpperBounds() {
    return this->mUpperBounds;
}

void DataConfigurations::CheckKernelValue(const string& aKernel) {

    // Check if the kernel name exists in the availableKernels set.
    if (availableKernels.count(aKernel) <= 0) {
        throw range_error("Invalid value for Kernel. Please check manual.");
    }
    else{
        // Check if the string is already in CamelCase format
        if (IsCamelCase(aKernel)) {
            this->SetKernel(aKernel);
            return ;
        }
        string str = aKernel;
        // Replace underscores with spaces and split the string into words
        std::replace(str.begin(), str.end(), '_', ' ');
        std::istringstream iss(str);
        std::string word, result;
        while (iss >> word) {
            // Capitalize the first letter of each word and append it to the result
            word[0] = toupper(word[0]);
            result += word;
        }
        this->SetKernel(result);
    }
}
bool DataConfigurations::IsCamelCase(std::string aString) {
    // If the string contains an underscore, it is not in CamelCase format
    if (aString.find('_') != std::string::npos) {
        return false;
    }
    // If the string starts with a lowercase letter, it is not in CamelCase format
    if (islower(aString[0])) {
        return false;
    }
    // If none of the above conditions hold, the string is in CamelCase format
    return true;
}

int DataConfigurations::ParseTheta(std::string aTheta) {
    if (aTheta.empty() || aTheta == ""){
        cout << "Empty?: " << this->GetParamtersNumber() << endl;
        for(int i = 0; i < this->GetParamtersNumber(); i++){
            std::cout << i << endl;
        }
    }
    return 0;
}
