
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file Results.cpp
 * @brief Defines the Results class for storing and accessing result data.
 * @version 1.0.1
 * @author Mahmoud ElKarargy
 * @date 2023-09-14
**/

#include <results/Results.hpp>
#include <common/Utils.hpp>

using namespace exageostat::results;
using namespace exageostat::common;

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

void Results::SetLoggerPath(const std::string &aLoggerPath) {
    this->mLoggerPath = aLoggerPath;
}

void Results::PrintEndSummary() {

    Verbose temp = exageostat::configurations::Configurations::GetVerbosity();
    exageostat::configurations::Configurations::SetVerbosity(STANDARD_MODE);
    LOGGER("")
    LOGGER("********************SUMMARY**********************")

    auto locations_number = this->mGeneratedLocationsNumber;
    if (locations_number > 0) {
        LOGGER("---- Data Generation Results ----")
        if (this->mIsSynthetic) {
            LOGGER(" #Synthetic Dataset")
        } else {
            LOGGER(" #Real Dataset")
        }
        LOGGER("  #Number of Locations: " << locations_number)
        if (this->mIsLogger && this->mIsSynthetic) {
            LOGGER("  #Data is written to file (", true)
            if (this->mLoggerPath.empty()) {
                this->mLoggerPath = LOG_PATH;
            }
            LOGGER_PRECISION(this->mLoggerPath << ").")
            LOGGER("")
        }
        LOGGER(" #Total Data Generation Execution Time: " << this->mExecutionTimeDataGeneration)
        LOGGER(" #Total Data Generation Gflop/s: " << this->mFlopsDataGeneration)
        LOGGER("")
    }
    if (this->mMLEIterations > 0) {
        LOGGER("---- Data Modeling Results ----")
        LOGGER(" #Number of MLE Iterations till reach Maximum: " << this->mMLEIterations)
        LOGGER(" #Found Maximum Theta at: ", true)
        for (double i: this->mMaximumTheta) {
            LOGGER_PRECISION(i << " ", 8)
        }
        LOGGER("")
        LOGGER(" #Final Log Likelihood value: " << this->mLogLikValue)
        LOGGER(" #Average Time Modeling per Iteration: " << this->GetAverageModelingExecutionTime())
        LOGGER(" #Average Flops per Iteration: " << this->GetAverageModelingFlops())
        LOGGER(" #Total MLE Execution time: " << this->mTotalModelingExecutionTime)
        LOGGER(" #Total MLE Gflop/s: " << this->mTotalModelingFlops)
        LOGGER("")
    }

    if (this->mZMiss > 0) {
        LOGGER("---- Data Prediction Results ----")
        LOGGER(" #Number of Missing Observations: " << this->mZMiss)
        if (this->mMSPEError > 0) {
            LOGGER(" #MSPE")
            LOGGER("  #MSPE Prediction Execution Time: " << this->mExecutionTimeMSPE)
            LOGGER("  #MSPE Gflop/s: " << this->mFlopsMSPE)
            LOGGER("  #Mean Square Error MSPE: " << this->mMSPEError)

        }
        if (!this->mIDWError.empty()) {
            LOGGER(" #IDW")
            LOGGER("  #IDW Error: ( ", true)
            for (int i = 0; i < 3; i++) {
                LOGGER_PRECISION(this->mIDWError[i] << " ", 8)
            }
            LOGGER_PRECISION(").")
            LOGGER("")
        }
        if (this->mMLOE > 0 || this->mMMOM > 0) {
            LOGGER(" #MLOE MMOM")
            LOGGER("  #MLOE: " << this->mMLOE)
            LOGGER("  #MMOM: " << this->mMMOM)
            LOGGER("  #MLOE-MMOM Execution Time: " << this->mExecutionTimeMLOEMMOM)
            LOGGER("  #MLOE-MMOM Matrix Generation Time: " << this->mGenerationTimeMLOEMMOM)
            LOGGER("  #MLOE-MMOM Cholesky Factorization Time: " << this->mFactoTimeMLOEMMOM)
            LOGGER("  #MLOE-MMOM Loop Time: " << this->mLoopTimeMLOEMMOM)
            LOGGER("  #MLOE-MMOM Number of flops: " << this->mFlopsMLOEMMOM)
            LOGGER("")
        }
    }
    if (this->mFisher00 != 0) {
        LOGGER(" #Fisher")
        LOGGER("  #Sd For Sigma2: " << this->mFisher00)
        LOGGER("  #Sd For Alpha: " << this->mFisher11)
        LOGGER("  #Sd For Nu: " << this->mFisher22)
        LOGGER("  #Fisher Execution Time: " << this->mTotalFisherTime)
        LOGGER("")
    }
    LOGGER("*************************************************")
    exageostat::configurations::Configurations::SetVerbosity(temp);
}

void Results::SetMLEIterations(int aIterationsNumber) {
    this->mMLEIterations = aIterationsNumber;
}

void Results::SetMaximumTheta(const std::vector<double> &aMaximumTheta) {
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

void Results::SetIDWError(const std::vector<double> &aIDWError) {
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
    throw std::runtime_error("Number of MLE Iterations is not set!");
}

double Results::GetAverageModelingFlops() const {
    if (this->mMLEIterations) {
        return this->mTotalModelingFlops / this->mMLEIterations;
    }
    throw std::runtime_error("Number of MLE Iterations is not set!");
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

void Results::SetFisher00(double aFisher00) {
    this->mFisher00 = aFisher00;
}

void Results::SetFisher11(double aFisher11) {
    this->mFisher11 = aFisher11;
}

void Results::SetFisher22(double aFisher22) {
    this->mFisher22 = aFisher22;
}


double Results::GetMLOE() const {
    return this->mMLOE;
}

double Results::GetMSPEError() const {
    return this->mMSPEError;
}

std::vector<double> Results::GetIDWError() const {
    return this->mIDWError;
}

double Results::GetMMOM() const {
    return this->mMMOM;
}

double Results::GetFisher00() const {
    return this->mFisher00;
}


double Results::GetFisher11() const {
    return this->mFisher11;
}

double Results::GetFisher22() const {
    return this->mFisher22;
}

