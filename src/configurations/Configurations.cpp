
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
using namespace exageostat::dataunits;
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

void Configurations::SetComputation(Computation aComputation) {
    this->mComputation = aComputation;
}

Computation Configurations::GetComputation() {
    return this->mComputation;
}

Computation Configurations::CheckComputationValue(std::string aValue) {

    if (aValue != "exact" and aValue != "Exact" and aValue != "Dense" and aValue != "dense"
        and aValue != "diag_approx" and aValue != "diagonal_approx"
        and aValue != "lr_approx" and aValue != "tlr" and aValue != "TLR") {
        throw range_error("Invalid value for Computation. Please use Exact, diagonal_approx or TLR.");
    }
    if (aValue == "exact" or aValue == "Exact" or aValue == "Dense" or aValue == "dense") {
        return EXACT_DENSE;
    } else if (aValue == "diag_approx" or aValue == "diagonal_approx") {
        return DIAGONAL_APPROX;
    }
    return TILE_LOW_RANK;
}
