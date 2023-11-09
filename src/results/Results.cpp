
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file Results.cpp
 * @brief Defines the Results class for storing and accessing result data.
 * @version 1.0.0
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

    auto locations_number = mGeneratedLocationsNumber;
    if (locations_number > 0) {
        LOGGER("---- Data Generation Results ----")
        if (mIsSynthetic) {
            LOGGER(" #Synthetic Dataset")
        } else {
            LOGGER(" #Real Dataset")
        }
        LOGGER("  #Number of Locations: " << locations_number)
        if (mIsLogger) {
            LOGGER("  #Data is written to file (", true)
            if (mLoggerPath.empty()) {
                mLoggerPath = LOG_PATH;
            }
            LOGGER_PRECISION(mLoggerPath << ").")
            LOGGER("")
        }
        LOGGER(" #Total Data Generation Execution Time: " << mExecutionTimeDataGeneration)
        LOGGER(" #Total Data Generation Gflop/s: " << mFlopsDataGeneration)
        LOGGER("")
    }
    if (mMLEIterations > 0) {
        LOGGER("---- Data Modeling Results ----")
        LOGGER(" #Number of MLE Iterations till reach Maximum: " << mMLEIterations)
        LOGGER(" #Found Maximum Theta at: ", true)
        for (double i: mMaximumTheta) {
            LOGGER_PRECISION(i << " ", 8)
        }
        LOGGER("")
        LOGGER(" #Final Log Likelihood value: " << mLogLikValue)
        LOGGER(" #Average Time Modeling per Iteration: " << this->GetAverageModelingExecutionTime())
        LOGGER(" #Average Flops per Iteration: " << this->GetAverageModelingFlops())
        LOGGER(" #Total MLE Execution time: " << mTotalModelingExecutionTime)
        LOGGER(" #Total MLE Gflop/s: " << mTotalModelingFlops)
        LOGGER("")
    }

    if (mZMiss > 0) {
        LOGGER("---- Data Prediction Results ----")
        LOGGER(" #Number of Missing Observations: " << mZMiss)
        if (mMSPEError > 0) {
            LOGGER(" #MSPE")
            LOGGER("  #MSPE Prediction Execution Time: " << mExecutionTimeMSPE)
            LOGGER("  #MSPE Gflop/s: " << mFlopsMSPE)
            LOGGER("  #Mean Square Error MSPE: " << mMSPEError)

        }
        if (!mIDWError.empty()) {
            LOGGER(" #IDW")
            LOGGER("  #IDW Error: ( ", true)
            for (int i = 0; i < 3; i++) {
                LOGGER_PRECISION(mIDWError[i] << " ", 8)
            }
            LOGGER_PRECISION(").")
            LOGGER("")
        }
        if (mMLOE > 0 || mMMOM > 0) {
            LOGGER(" #MLOE MMOM")
            LOGGER("  #MLOE: " << mMLOE)
            LOGGER("  #MMOM: " << mMMOM)
            LOGGER("  #MLOE-MMOM Execution Time: " << mExecutionTimeMLOEMMOM)
            LOGGER("  #MLOE-MMOM Matrix Generation Time: " << mGenerationTimeMLOEMMOM)
            LOGGER("  #MLOE-MMOM Cholesky Factorization Time: " << mFactoTimeMLOEMMOM)
            LOGGER("  #MLOE-MMOM Loop Time: " << mLoopTimeMLOEMMOM)
            LOGGER("  #MLOE-MMOM Number of flops: " << mFlopsMLOEMMOM)
            LOGGER("")
        }
    }
    if(mFisher00 != 0){
        LOGGER(" #Fisher")
        LOGGER("  #Sd For Sigma2: " << mFisher00)
        LOGGER("  #Sd For Alpha: " << mFisher11)
        LOGGER("  #Sd For Nu: " << mFisher22)
        LOGGER("  #Fisher Execution Time: " << mTotalFisherTime)
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