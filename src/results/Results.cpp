
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
        }
        LOGGER("")
    }
    if (mMLEIterations > 0) {
        LOGGER("---- Data Modeling Results ----")
        LOGGER(" #Number of MLE Iterations till reach Maximum: " << mMLEIterations)
        LOGGER(" #Found Maximum Theta at ( ", true)
        for (double i: mMaximumTheta) {
            LOGGER_PRECISION(i << " ", 8)
        }
        LOGGER_PRECISION(").")
        LOGGER("")
        LOGGER(" #Maximum Log Likelihood value: " << mLogLikValue)
    }
    if (mZMiss > 0) {
        LOGGER("---- Data Prediction Results ----")
        LOGGER(" #Number of Missing Observations: " << mZMiss)
        if (mMSPEError > 0) {
            LOGGER(" #MSPE")
            LOGGER("  # Prediction Execution Time: " << mExecutionTime)
            LOGGER("  # Flops: " << mFlops)
            LOGGER("  # Mean Square Error (MSE): " << mMSPEError)

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
        }
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

void Results::SetExecutionTime(double aExecutionTime) {
    this->mExecutionTime = aExecutionTime;
}

void Results::SetFlops(double aFlops) {
    this->mFlops = aFlops;
}

Results *Results::mpInstance = nullptr;
