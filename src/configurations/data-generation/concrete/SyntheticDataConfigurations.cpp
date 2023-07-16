
/*
 * Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
 * All rights reserved.
 * ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).
 */

/**
 * @file SyntheticDataConfigurations.cpp
 * @brief Implementation for Synthetic data Configurations.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-02-01
**/

#include <iostream>

#include <configurations/data-generation/concrete/SyntheticDataConfigurations.hpp>

using namespace std;

using namespace exageostat::configurations::data_configurations;
using namespace exageostat::common;


void SyntheticDataConfigurations::InitModuleArguments(int aArgC, char **apArgV) {

    InitializeArguments(aArgC, apArgV);

    string argument;
    string argument_name;
    string argument_value;
    int equal_sign_Idx;

    // Loop through the arguments that are specific for data generation.
    for (int i = 1; i < aArgC; ++i) {
        argument = apArgV[i];
        equal_sign_Idx = static_cast<int>(argument.find('='));
        argument_name = argument.substr(0, equal_sign_Idx);

        // Check if argument has an equal sign.
        if (equal_sign_Idx != string::npos) {
            argument_value = argument.substr(equal_sign_Idx + 1);

            // Check the argument name and set the corresponding value
            if (argument_name == "--Dimension" || argument_name == "--dimension"
                       || argument_name == "--dim" || argument_name == "--Dim") {
                SetDimension(CheckDimensionValue(argument_value));
            } else if (argument_name == "--ZmissNumber" || argument_name == "--Zmiss") {
                SetUnknownObservationsNb(CheckUnknownObservationsValue(argument_value));
            } else if (argument_name == "--lb" || argument_name == "--olb" || argument_name == "--lower_bounds") {
                std::vector<double> theta = ParseTheta(argument_value);
                SetLowerBounds(theta);
            } else if (argument_name == "--ub" || argument_name == "--oub" || argument_name == "--upper_bounds") {
                std::vector<double> theta = ParseTheta(argument_value);
                SetUpperBounds(theta);
            } else if (argument_name == "--target_theta" || argument_name == "--ttheta" ||
                       argument_name == "--tTheta") {
                std::vector<double> theta = ParseTheta(argument_value);
                SetTargetTheta(theta);
            }
        } else {
            if (argument_name == "--syntheticData" || argument_name == "--SyntheticData" ||
                argument_name == "--synthetic_data" || argument_name == "--synthetic") {
                SetIsSynthetic(true);
            }
        }
    }
    // Throw Errors if any of these arguments aren't given by the user.
    if (GetKernel().empty()) {
        throw domain_error("You need to set the Kernel, before starting");
    }
}

SyntheticDataConfigurations::SyntheticDataConfigurations(int aArgC, char **apArgV) {

    InitializeArguments(aArgC, apArgV);
    InitModuleArguments(aArgC, apArgV);
}

Dimension SyntheticDataConfigurations::GetDimension() {
    return this->mDimension;
}

void SyntheticDataConfigurations::SetDimension(Dimension aDimension) {
    this->mDimension = aDimension;
}

Dimension SyntheticDataConfigurations::CheckDimensionValue(const string &aDimension) {

    if (aDimension != "2D" and aDimension != "2d"
        and aDimension != "3D" and aDimension != "3d"
        and aDimension != "st" and aDimension != "ST") {
        throw range_error("Invalid value for Dimension. Please use 2D, 3D or ST.");
    }
    if (aDimension == "2D" or aDimension == "2d") {
        return Dimension2D;
    } else if (aDimension == "3D" or aDimension == "3d") {
        return Dimension3D;
    }
    return DimensionST;
}

int SyntheticDataConfigurations::CheckUnknownObservationsValue(const string &aValue) {
    int value = CheckNumericalValue(aValue);
    if (value >= GetProblemSize()) {
        throw range_error("Invalid value for ZmissNumber. Please make sure it's smaller than Problem size");
    }
    return value;
}

