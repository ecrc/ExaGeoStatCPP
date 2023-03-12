
/*
 * Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
 * Copyright (C) 2023 by Brightskies inc,
 * All rights reserved.
 * ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).
 */

/**
 * @file Configurations.cpp
 *
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-01-31
**/

#include <iostream>
#include <configurations/Configurations.hpp>

using namespace exageostat::configurations;
using namespace std;


int Configurations::GetProblemSize() {
    return this->mProblemSize;
}

void Configurations::SetProblemSize(int aProblemSize) {
    this->mProblemSize = aProblemSize;
}

int Configurations::CheckNumericalValue(string aValue) {

    int numericalValue = -1;
    try {
        numericalValue = stoi(aValue);
    }
    catch (...) {
        throw range_error("Invalid value. Please use Numerical values only.");
    }

    if (numericalValue < 0) {
        throw range_error("Invalid value. Please use positive values");
    }
    return numericalValue;
}

void Configurations::SetTimeSlot(int aTimeSlot) {
    this->mTimeSlot = aTimeSlot;
}

int Configurations::GetTimeSlot() {
    return this->mTimeSlot;
}
