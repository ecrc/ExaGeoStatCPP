
/*
 * Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
 * Copyright (C) 2023 by Brightskies inc,
 * All rights reserved.
 * ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).
 */

/**
 * @file Configurations.cpp
 * @brief This file defines the Configurations class which stores the configuration parameters for ExaGeoStat.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-01-31
**/

#include <iostream>
#include <configurations/Configurations.hpp>
#include <utility>

using namespace exageostat::configurations;
using namespace exageostat::common;
using namespace std;


int Configurations::GetProblemSize() const {
    return this->mProblemSize;
}

void Configurations::SetProblemSize(int aProblemSize) {
    this->mProblemSize = aProblemSize;
}

void Configurations::SetTimeSlot(int aTimeSlot) {
    this->mTimeSlot = aTimeSlot;
}

int Configurations::GetTimeSlot() const {
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

void Configurations::SetOperator(Operators aOperator) {
    this->mOperator = aOperator;
}

Operators Configurations::GetOperator() {
    return this->mOperator;
}

void Configurations::SetPGrid(int aPGrid) {
    this->mPGrid = aPGrid;
}

int Configurations::GetPGrid() const {
    return this->mPGrid;
}

void Configurations::SetP(int aP) {
    this->mP = aP;
}

int Configurations::GetP() const {
    return this->mP;
}

void Configurations::SetDenseTileSize(int aTileSize) {
    this->mDenseTileSize = aTileSize;
}

int Configurations::GetDenseTileSize() const {
    return this->mDenseTileSize;
}

void Configurations::SetLowTileSize(int aTileSize) {
    this->mLowTileSize = aTileSize;
}

int Configurations::GetLowTileSize() const {
    return this->mLowTileSize;
}

void Configurations::SetQGrid(int aQGrid) {
    this->mQGrid = aQGrid;
}

int Configurations::GetQGrid() const {
    return this->mQGrid;

}

std::vector<void *> &Configurations::GetDescriptorC() {
    return this->mpDescriptorC;
}

std::vector<void *> &Configurations::GetDescriptorZ() {
    return this->mpDescriptorZ;
}

void *&Configurations::GetDescriptorZcpy() {
    return this->mpDescriptorZcpy;
}

std::vector<void *> &Configurations::GetDescriptorProduct() {
    return this->mpDescriptorProduct;
}

void *&Configurations::GetDescriptorDeterminant() {
    return this->mpDescriptorDeterminant;
}

std::vector<void *> &Configurations::GetDescriptorCD() {
    return this->mpDescriptorCD;
}

std::vector<void *> &Configurations::GetDescriptorCUV() {
    return this->mpDescriptorCUV;
}

std::vector<void *> &Configurations::GetDescriptorCrk() {
    return this->mpDescriptorCrk;
}

void *&Configurations::GetDescriptorZObservations() {
    return this->mpDescriptorZObservations;
}

void *&Configurations::GetDescriptorMSE() {
    return this->mpDescriptorMSE;
}

void *&Configurations::GetDescriptorZActual() {
    return this->mpDescriptorZActual;
}

void Configurations::SetCoresNumber(int aCoresNumbers) {
    this->mCoresNumber = aCoresNumbers;
}

int Configurations::GetCoresNumber() const {
    return this->mCoresNumber;
}

void Configurations::SetGPUsNumber(int aGPUsNumber) {
    this->mGPUsNumber = aGPUsNumber;
}

int Configurations::GetGPUsNumber() const {
    return this->mGPUsNumber;
}

void Configurations::SetIsOOC(bool aIsOOC) {
    this->mIsOOC = aIsOOC;
}

bool Configurations::GetIsOOC() const {
    return this->mIsOOC;
}

void Configurations::SetMaxRank(int aMaxRank) {
    this->mMaxRank = aMaxRank;
}

int Configurations::GetMaxRank() {
    return this->mMaxRank;
}

void Configurations::SetUnknownObservationsNb(int aUnknownObservationsNumber) {
    this->mUnknownObservationsNumber = aUnknownObservationsNumber;
}

int Configurations::GetUnknownObservationsNb() const {
    return this->mUnknownObservationsNumber;
}

void Configurations::SetKnownObservationsValues(int aKnownObservationsValues) {
    this->mKnownObservationsValues = aKnownObservationsValues;
}

int Configurations::GetKnownObservationsValues() {
    return this->mKnownObservationsValues;
}

int Configurations::GetApproximationMode() const {
    return this->mApproximationMode;
}

void Configurations::SetApproximationMode(int aApproximationMode) {
    this->mApproximationMode = aApproximationMode;
}

double Configurations::GetMeanSquareError() {
    return this->mMeanSquareError;
}

void Configurations::SetMeanSquareError(double aMeanSquareError) {
    this->mMeanSquareError = aMeanSquareError;
}

void Configurations::SetActualObservationsFilePath(std::string aKnownObservationsValues) {
    this->mActualObservationsFilePath = std::move(aKnownObservationsValues);
}

string Configurations::GetActualObservationsFilePath() {
    return this->mActualObservationsFilePath;
}

void Configurations::SetDeterminantValue(double aDeterminantValue) {
    this->mDeterminantValue = aDeterminantValue;
}

double Configurations::GetDeterminantValue() {
    return this->mDeterminantValue;
}

int Configurations::GetSeed() {
    return this->mSeed;
}

void Configurations::SetSeed(int aSeed) {
    this->mSeed = aSeed;
}

int Configurations::CheckNumericalValue(const string &aValue) {

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

Computation Configurations::CheckComputationValue(const std::string &aValue) {

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

Precision Configurations::CheckPrecisionValue(const std::string &aValue) {

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

void Configurations::SetSequence(void *apSequence) {
    this->mpSequence = apSequence;
}

void *Configurations::GetSequence() {
    return this->mpSequence;
}

void Configurations::SetRequest(void *apRequest) {
    this->mpRequest = apRequest;
}

void *Configurations::GetRequest() {
    return this->mpRequest;
}

RunMode Configurations::mRunMode = RunMode::STANDARD_MODE;

RunMode Configurations::GetRunMode() {
    return Configurations::mRunMode;
}

void Configurations::SetRunMode(RunMode aRunMode) {
    Configurations::mRunMode = aRunMode;
}

void Configurations::SetLogger(bool aLogger) {
    this->mLogger = aLogger;
}

bool Configurations::GetLogger() {
    return this->mLogger;
}

std::string *Configurations::GetLoggerPath() {
    return &this->mLoggerPath;
}

void Configurations::SetLoggerPath(const string &aLoggerPath) {
    this->mLoggerPath = aLoggerPath;
}
