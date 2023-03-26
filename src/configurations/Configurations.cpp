
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
using namespace exageostat::common;
using namespace std;


int Configurations::GetProblemSize() {
    return this->mProblemSize;
}

void Configurations::SetProblemSize(int aProblemSize) {
    this->mProblemSize = aProblemSize;
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

Precision Configurations::GetPrecision() {
    return this->mPrecision;
}

void Configurations::SetPrecision(Precision aPrecision) {
    this->mPrecision = aPrecision;
}

void Configurations::SetPGrid(int aPGrid) {
    this->mPGrid = aPGrid;
}

int Configurations::GetPGrid() {
    return this->mPGrid;
}

void Configurations::SetP(int aP) {
    this->mP = aP;
}

int Configurations::GetP() {
    return this->mP;
}

void Configurations::SetDenseTileSize(int aTileSize) {
    this->mDenseTileSize = aTileSize;
}

int Configurations::GetDenseTileSize() {
    return this->mDenseTileSize;
}

void Configurations::SetLowTileSize(int aTileSize) {
    this->mLowTileSize = aTileSize;
}

int Configurations::GetLowTileSize() {
    return this->mLowTileSize;
}

void Configurations::SetQGrid(int aQGrid) {
    this->mQGrid = aQGrid;
}

int Configurations::GetQGrid() {
    return this->mQGrid;

}
std::vector<void *> Configurations::GetDescriptorC() {
    return this->mpDescriptorC;
}
std::vector<void *> Configurations::GetDescriptorZ() {
    return this->mpDescriptorZ;
}
void * Configurations::GetDescriptorZcpy() {
    return this->mpDescriptorZcpy;
}
std::vector<void *> Configurations::GetDescriptorProduct() {
    return this->mpDescriptorProduct;
}
void * Configurations::GetDescriptorDeterminant() {
    return this->mpDescriptorDeterminant;
}

void Configurations::SetCoresNumber(int aCoresNumbers) {
    this->mCoresNumber = aCoresNumbers;
}

int *Configurations::GetCoresNumber() {
    return &(this->mCoresNumber);
}

void Configurations::SetGPUsNumber(int aGPUsNumber) {
    this->mGPUsNumber = aGPUsNumber;
}

int *Configurations::GetGPUsNumber() {
    return &(this->mGPUsNumber);
}

void Configurations::SetIsOOC(int aIsOOC) {
    this->mIsOOC = aIsOOC;
}

int Configurations::GetIsOOC() {
    return this->mIsOOC;
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

Precision Configurations::CheckPrecisionValue(std::string aValue) {

    if (aValue != "single" and aValue != "Single"
        and aValue != "double" and aValue != "Double"
        and aValue != "mix" and aValue != "Mix" and aValue != "Mixed" and aValue != "mixed") {
        throw range_error("Invalid value for Computation. Please use Single, Double or Mixed.");
    }
    if (aValue == "single" or aValue == "Single") {
        return SINGLE;
    } else if (aValue == "double" or aValue == "Double") {
        return DOUBLE;
    }
    return MIXED;

}
