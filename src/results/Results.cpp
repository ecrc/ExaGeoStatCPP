
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
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

Results *Results::mpInstance = nullptr;

Results *Results::GetInstance() {

    if (mpInstance == nullptr) {
        mpInstance = new Results();
    }
    return mpInstance;
}

void Results::UpdateDictionary(const string &key, const string &value) {
    mSummaryDictionary[key] = value;
}

void Results::SetIsSynthetic(bool aIsSynthetic, const std::string &aKey) {
    this->mIsSynthetic = aIsSynthetic;
    UpdateDictionary(aKey.empty() ? "Data is Synthetic" : aKey, aIsSynthetic ? "Yes" : "No");
}

void Results::SetGeneratedLocationsNumber(int aNumLocations, const std::string &aKey) {
    this->mGeneratedLocationsNumber = aNumLocations;
    UpdateDictionary(aKey.empty() ? "Number of Locations" : aKey, to_string(aNumLocations));
}

void Results::SetIsLogger(bool aIsLogger, const std::string &aKey) {
    this->mIsLogger = aIsLogger;
    UpdateDictionary(aKey.empty() ? "Logger Enabled" : aKey, aIsLogger ? "Yes" : "No");
}

void Results::SetLoggerPath(const string &aLoggerPath, const std::string &aKey) {
    this->mLoggerPath = aLoggerPath;
    UpdateDictionary(aKey.empty() ? "Logger Path" : aKey, aLoggerPath.empty() ? LOG_PATH : aLoggerPath);
}

void Results::PrintEndSummary() {

    Verbose temp = Configurations::GetVerbosity();
    Configurations::SetVerbosity(STANDARD_MODE);
    LOGGER("********************SUMMARY**********************")
    for (const auto &entry: mSummaryDictionary) {
        LOGGER("#" << entry.first << ": " << entry.second)
    }
    LOGGER("*************************************************")
    Configurations::SetVerbosity(temp);
}

void Results::SetMLEIterations(int aIterationsNumber, const std::string& aKey) {
    this->mMLEIterations = aIterationsNumber;
    UpdateDictionary(aKey.empty() ? "Number of MLE Iterations" : aKey, to_string(aIterationsNumber));
}

void Results::SetMaximumTheta(const vector<double> &aMaximumTheta, const std::string& aKey) {
    this->mMaximumTheta = aMaximumTheta;

    ostringstream oss;
    oss << "[ ";
    for (double val: aMaximumTheta) {
        oss << fixed << setprecision(8) << val << " ";
    }
    oss << "]";
    UpdateDictionary(aKey.empty() ? "Found Maximum Theta at" : aKey, oss.str());
}

void Results::SetLogLikValue(double aLogLikValue, const std::string& aKey) {
    this->mLogLikValue = aLogLikValue;
    UpdateDictionary(aKey.empty() ? "Final Log Likelihood Value" : aKey, to_string(aLogLikValue));
}

void Results::SetZMiss(int aZMiss, const std::string& aKey) {
    this->mZMiss = aZMiss;
    UpdateDictionary(aKey.empty() ? "Number of Missing Observations" : aKey, to_string(aZMiss));
}

void Results::SetMSPEError(double aMSPEError, const std::string& aKey) {
    this->mMSPEError = aMSPEError;
    UpdateDictionary(aKey.empty() ? "Mean Square Error MSPE" : aKey, to_string(aMSPEError));
}

void Results::SetIDWError(const vector<double> &aIDWError, const std::string& aKey) {
    this->mIDWError = aIDWError;

    ostringstream oss;
    oss << "[ ";
    for (double val: aIDWError) {
        oss << fixed << setprecision(8) << val << " ";
    }
    oss << "]";
    UpdateDictionary(aKey.empty() ? "IDW Error" : aKey, oss.str());
}

void Results::SetMLOE(double aMLOE, const std::string& aKey) {
    this->mMLOE = aMLOE;
    UpdateDictionary(aKey.empty() ? "MLOE" : aKey, to_string(aMLOE));
}

void Results::SetMMOM(double aMMOM, const std::string& aKey) {
    this->mMMOM = aMMOM;
    UpdateDictionary(aKey.empty() ? "MMOM" : aKey, to_string(aMMOM));
}

void Results::SetMSPEExecutionTime(double aTime, const std::string& aKey) {
    this->mExecutionTimeMSPE = aTime;
    UpdateDictionary(aKey.empty() ? "MSPE Prediction Execution Time" : aKey, to_string(aTime));
}

void Results::SetMSPEFlops(double aFlops, const std::string& aKey) {
    this->mFlopsMSPE = aFlops;
    UpdateDictionary(aKey.empty() ? "MSPE Gflop/s" : aKey, to_string(aFlops));
}

void Results::SetTotalModelingExecutionTime(double aTime, const std::string& aKey) {
    this->mTotalModelingExecutionTime = aTime;
    UpdateDictionary(aKey.empty() ? "Total Modeling Execution Time" : aKey, to_string(aTime));
    if (this->mMLEIterations) {
        UpdateDictionary("Average Time Modeling per Iteration", to_string(
        this->mTotalModelingExecutionTime / this->mMLEIterations));
    }
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

void Results::SetTotalModelingFlops(double aTime, const std::string& aKey) {
    this->mTotalModelingFlops = aTime;
    UpdateDictionary(aKey.empty() ? "Total Modeling Flops" : aKey, to_string(aTime));
    if (this->mMLEIterations) {
        UpdateDictionary("Average Flops per Iteration", to_string(this->mTotalModelingFlops / this->mMLEIterations));
    }
}

double Results::GetTotalModelingFlops() const {
    return this->mTotalModelingFlops;
}

void Results::SetExecutionTimeMLOEMMOM(double aTime, const std::string& aKey) {
    this->mExecutionTimeMLOEMMOM = aTime;
    UpdateDictionary(aKey.empty() ? "MLOE-MMOM Execution Time" : aKey, to_string(aTime));
}

void Results::SetMatrixGenerationTimeMLOEMMOM(double aTime, const std::string& aKey) {
    this->mGenerationTimeMLOEMMOM = aTime;
    UpdateDictionary(aKey.empty() ? "MLOE-MMOM Matrix Generation Time" : aKey, to_string(aTime));

}

void Results::SetFactoTimeMLOEMMOM(double aTime, const std::string& aKey) {
    this->mFactoTimeMLOEMMOM = aTime;
    UpdateDictionary(aKey.empty() ? "MLOE-MMOM Cholesky Factorization Time" : aKey, to_string(aTime));
}

void Results::SetLoopTimeMLOEMMOM(double aTime, const std::string& aKey) {
    this->mLoopTimeMLOEMMOM = aTime;
    UpdateDictionary(aKey.empty() ? "MLOE-MMOM Loop Time" : aKey, to_string(aTime));
}

void Results::SetFlopsMLOEMMOM(double aFlops, const std::string& aKey) {
    this->mFlopsMLOEMMOM = aFlops;
    UpdateDictionary(aKey.empty() ? "MLOE-MMOM Number of flops" : aKey, to_string(aFlops));
}

void Results::SetTotalDataGenerationExecutionTime(double aTime, const std::string& aKey) {
    this->mExecutionTimeDataGeneration = aTime;
    UpdateDictionary(aKey.empty() ? "Total Data Generation Execution Time" : aKey, to_string(aTime));
}

void Results::SetTotalDataGenerationFlops(double aFlops, const std::string& aKey) {
    this->mFlopsDataGeneration = aFlops;
    UpdateDictionary(aKey.empty() ? "Total Data Generation Gflop/s" : aKey, to_string(aFlops));
}

void Results::SetTotalFisherTime(double aTime, const std::string& aKey) {
    this->mTotalFisherTime = aTime;
    UpdateDictionary(aKey.empty() ? "Fisher Execution Time" : aKey, to_string(aTime));
}

void Results::SetFisherMatrix(vector<double> aFisherMatrix, const std::string& aKey) {
    this->mFisherMatrix = std::move(aFisherMatrix);

    ostringstream oss;
    if (this->mFisherMatrix.size() >= 3) {
        oss << "Sd For Sigma2: " << this->mFisherMatrix[0] << ", "
            << "Sd For Alpha: " << this->mFisherMatrix[1] << ", "
            << "Sd For Nu: " << this->mFisherMatrix[2];
    }
    UpdateDictionary(aKey.empty() ? "Fisher Matrix" : aKey, oss.str());
}

void Results::SetPredictedMissedValues(vector<double> aPredictedValues, const std::string& aKey) {
    this->mPredictedMissedValues = std::move(aPredictedValues);

    ostringstream oss;
    if (this->mPredictedMissedValues.size() >= 3) {
        oss << this->mPredictedMissedValues[0] << ", "
            << this->mPredictedMissedValues[1] << ", "
            << this->mPredictedMissedValues[2];
    }
    UpdateDictionary(aKey.empty() ? "Predicted Values" : aKey, oss.str());
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

Results *Results::mpInstance = nullptr;
