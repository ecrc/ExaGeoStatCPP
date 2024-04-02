
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file Results.cpp
 * @brief Defines the Results class for storing and accessing result data.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @date 2024-02-04
**/

#include <results/Results.hpp>
#include <utilities/Logger.hpp>
#include <utility>

using namespace exageostat::results;
using namespace exageostat::common;
using namespace exageostat::configurations;

using namespace std;

Results *Results::GetInstance() {

    if (mpInstance == nullptr) {
        mpInstance = new Results();
    }
    return mpInstance;
}

void Results::SetIsSynthetic(bool aIsSynthetic) {
    this->mIsSynthetic = aIsSynthetic;
}

void Results::SetGeneratedLocationsNumber(int aNumLocations) {
    this->mGeneratedLocationsNumber = aNumLocations;
}

void Results::SetIsLogger(bool aIsLogger) {
    this->mIsLogger = aIsLogger;
}

void Results::SetLoggerPath(const string &aLoggerPath) {
    this->mLoggerPath = aLoggerPath;
}

void Results::PrintEndSummary() {

    Verbose temp = Configurations::GetVerbosity();
    Configurations::SetVerbosity(STANDARD_MODE);
    LOGGER("********************SUMMARY**********************")

    auto locations_number = this->mGeneratedLocationsNumber;
    if (locations_number > 0) {
        LOGGER("#Number of Locations: " << locations_number)
        if (this->mIsLogger && this->mIsSynthetic) {
            LOGGER("  #Data is written to file (", true)
            if (this->mLoggerPath.empty()) {
                this->mLoggerPath = LOG_PATH;
            }
            LOGGER_PRECISION(this->mLoggerPath << ").")
            LOGGER("")
        }
        VERBOSE("#Total Data Generation Execution Time: " << this->mExecutionTimeDataGeneration)
        VERBOSE("#Total Data Generation Gflop/s: " << this->mFlopsDataGeneration)
    }
    if (this->mMLEIterations > 0) {
        LOGGER("#Number of MLE Iterations: " << this->mMLEIterations)
        LOGGER("#Found Maximum Theta at: ", true)
        for (double i: this->mMaximumTheta) {
            LOGGER_PRECISION(i << " ", 8)
        }
        LOGGER("")
        LOGGER("#Final Log Likelihood value: " << this->mLogLikValue)
        VERBOSE("#Average Time Modeling per Iteration: " << this->GetAverageModelingExecutionTime())
        VERBOSE("#Average Flops per Iteration: " << this->GetAverageModelingFlops())
        VERBOSE("#Total MLE Execution time: " << this->mTotalModelingExecutionTime)
        VERBOSE("#Total MLE GFlop/s: " << this->mTotalModelingFlops)
    }

    if (this->mZMiss > 0) {
        LOGGER("#Number of Missing Observations: " << this->mZMiss)
        if (this->mMSPEError > 0) {
            VERBOSE("#MSPE Prediction Execution Time: " << this->mExecutionTimeMSPE)
            VERBOSE("#MSPE Gflop/s: " << this->mFlopsMSPE)
            LOGGER("#Mean Square Error MSPE: " << this->mMSPEError)

        }
        if (!this->mIDWError.empty()) {
            LOGGER("#IDW Error: ( ", true)
            for (int i = 0; i < 3; i++) {
                LOGGER_PRECISION(this->mIDWError[i] << " ", 8)
            }
            LOGGER_PRECISION(").")
            LOGGER("")
        }
        if (this->mMLOE > 0 || this->mMMOM > 0) {
            LOGGER("#MLOE: " << this->mMLOE << "\t\t#MMOM: " << this->mMMOM)
            VERBOSE("#MLOE-MMOM Execution Time: " << this->mExecutionTimeMLOEMMOM)
            VERBOSE("#MLOE-MMOM Matrix Generation Time: " << this->mGenerationTimeMLOEMMOM)
            VERBOSE("#MLOE-MMOM Cholesky Factorization Time: " << this->mFactoTimeMLOEMMOM)
            VERBOSE("#MLOE-MMOM Loop Time: " << this->mLoopTimeMLOEMMOM)
            VERBOSE("#MLOE-MMOM Number of flops: " << this->mFlopsMLOEMMOM)
        }
    }
    if (!this->mFisherMatrix.empty()) {
        LOGGER("#Sd For Sigma2: " << this->mFisherMatrix[0])
        LOGGER("#Sd For Alpha: " << this->mFisherMatrix[1])
        LOGGER("#Sd For Nu: " << this->mFisherMatrix[2])
        VERBOSE("#Fisher Execution Time: " << this->mTotalFisherTime)
    }
    LOGGER("*************************************************")
    Configurations::SetVerbosity(temp);
}

void Results::SetMLEIterations(int aIterationsNumber) {
    this->mMLEIterations = aIterationsNumber;
}

void Results::SetMaximumTheta(const vector<double> &aMaximumTheta) {
    this->mMaximumTheta = aMaximumTheta;
}

void Results::SetLogLikValue(double aLogLikValue) {
    this->mLogLikValue = aLogLikValue;
}

void Results::SetZMiss(int aZMiss) {
    this->mZMiss = aZMiss;
}

void Results::SetMSPEError(double aMSPEError) {
    this->mMSPEError = aMSPEError;
}

void Results::SetIDWError(const vector<double> &aIDWError) {
    this->mIDWError = aIDWError;
}

void Results::SetMLOE(double aMLOE) {
    this->mMLOE = aMLOE;
}

void Results::SetMMOM(double aMMOM) {
    this->mMMOM = aMMOM;
}

void Results::SetMSPEExecutionTime(double aTime) {
    this->mExecutionTimeMSPE = aTime;
}

void Results::SetMSPEFlops(double aFlops) {
    this->mFlopsMSPE = aFlops;
}

void Results::SetTotalModelingExecutionTime(double aTime) {
    this->mTotalModelingExecutionTime = aTime;
}

double Results::GetTotalModelingExecutionTime() const {
    return this->mTotalModelingExecutionTime;
}

double Results::GetAverageModelingExecutionTime() const {
    if (this->mMLEIterations) {
        return this->mTotalModelingExecutionTime / this->mMLEIterations;
    }
    throw runtime_error("Number of MLE Iterations is not set!");
}

double Results::GetAverageModelingFlops() const {
    if (this->mMLEIterations) {
        return this->mTotalModelingFlops / this->mMLEIterations;
    }
    throw runtime_error("Number of MLE Iterations is not set!");
}

void Results::SetTotalModelingFlops(double aTime) {
    this->mTotalModelingFlops = aTime;
}

double Results::GetTotalModelingFlops() const {
    return this->mTotalModelingFlops;
}

Results *Results::mpInstance = nullptr;

void Results::SetExecutionTimeMLOEMMOM(double aTime) {
    this->mExecutionTimeMLOEMMOM = aTime;
}

void Results::SetMatrixGenerationTimeMLOEMMOM(double aTime) {
    this->mGenerationTimeMLOEMMOM = aTime;
}

void Results::SetFactoTimeMLOEMMOM(double aTime) {
    this->mFactoTimeMLOEMMOM = aTime;
}

void Results::SetLoopTimeMLOEMMOM(double aTime) {
    this->mLoopTimeMLOEMMOM = aTime;
}

void Results::SetFlopsMLOEMMOM(double aFlops) {
    this->mFlopsMLOEMMOM = aFlops;
}

void Results::SetTotalDataGenerationExecutionTime(double aTime) {
    this->mExecutionTimeDataGeneration = aTime;
}

void Results::SetTotalDataGenerationFlops(double aFlops) {
    this->mFlopsDataGeneration = aFlops;
}

void Results::SetTotalFisherTime(double aTime) {
    this->mTotalFisherTime = aTime;
}

void Results::SetFisherMatrix(vector<double> aFisherMatrix) {
    this->mFisherMatrix = std::move(aFisherMatrix);
}

void Results::SetPredictedMissedValues(vector<double> aPredictedValues) {
    this->mPredictedMissedValues = std::move(aPredictedValues);
}

double Results::GetMLOE() const {
    return this->mMLOE;
}

double Results::GetMSPEError() const {
    return this->mMSPEError;
}

vector<double> Results::GetIDWError() const {
    return this->mIDWError;
}

double Results::GetMMOM() const {
    return this->mMMOM;
}

std::vector<double> Results::GetFisherMatrix() const {
    return this->mFisherMatrix;
}

std::vector<double> Results::GetPredictedMissedValues() const {
    return this->mPredictedMissedValues;
}