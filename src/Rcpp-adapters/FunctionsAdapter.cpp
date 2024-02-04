
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file FunctionsAdapter.cpp
 * @brief Header file for function adapters in the ExaGeoStat software.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @date 2024-01-29
**/

#include <cstring>

#include <Rcpp-adapters/FunctionsAdapter.hpp>
#include <data-units/ExaGeoStatData.hpp>
#include <api/ExaGeoStat.hpp>

using namespace std;

using namespace exageostat::dataunits;
using namespace exageostat::api;

namespace exageostat::adapters {

    Configurations *
    R_InitializeArguments(const int &aProblemSize, const string &aKernelName, const vector<int> &aTileSize,
                         const vector<int> &aP_QGrid, const int &aTimeSlot, const string &aComputation,
                         const string &aPrecision, const vector<int> &aCoresGPUsNumber, const int &aBand,
                         const int &aMaxRank, const vector<double> &aInitialTheta,
                         const vector <vector<double>> &aLowerUpperBounds, const vector<double> &aEstimatedTheta,
                         const string &aVerbose, const string &aDimension, const int &aMaxMleIterations,
                         const double &aTolerance, const vector<int> &aPrediction) {

        vector<string> argStrings;

        // Program name or placeholder
        argStrings.emplace_back("--ConfigurationsModuleR");
        // Convert each argument to string and add to argv
        argStrings.push_back("--N=" + to_string(aProblemSize));
        argStrings.push_back("--kernel=" + aKernelName);
        argStrings.push_back("--dts=" + to_string(aTileSize[0]));
        argStrings.push_back("--lts=" + to_string(aTileSize[1]));
        argStrings.push_back("--p=" + to_string(aP_QGrid[0]));
        argStrings.push_back("--q=" + to_string(aP_QGrid[1]));
        argStrings.push_back("--time_slot=" + to_string(aTimeSlot));
        argStrings.push_back("--computation=" + aComputation);
        argStrings.push_back("--precision=" + aPrecision);
        argStrings.push_back("--cores=" + to_string(aCoresGPUsNumber[0]));
        argStrings.push_back("--gpus=" + to_string(aCoresGPUsNumber[1]));
        argStrings.push_back("--band=" + to_string(aBand));
        argStrings.push_back("--max_rank=" + to_string(aMaxRank));
        argStrings.push_back("--iTheta=" + to_string(aInitialTheta[0]) + ":" + to_string(aInitialTheta[1]) + ":" +
                       to_string(aInitialTheta[2]));
        argStrings.push_back("--lb=" + to_string(aLowerUpperBounds[0][0]) + ":" + to_string(aLowerUpperBounds[0][1]) + ":" +
                       to_string(aLowerUpperBounds[0][2]));
        argStrings.push_back("--ub=" + to_string(aLowerUpperBounds[1][0]) + ":" + to_string(aLowerUpperBounds[1][1]) + ":" +
                       to_string(aLowerUpperBounds[1][2]));
        argStrings.push_back("--eTheta=" + to_string(aEstimatedTheta[0]) + ":" + to_string(aEstimatedTheta[1]) + ":" +
                       to_string(aEstimatedTheta[2]));
        argStrings.push_back("--verbose=" + aVerbose);
        argStrings.push_back("--dimension=" + aDimension);
        argStrings.push_back("--max_mle_iterations=" + to_string(aMaxMleIterations));
        argStrings.push_back("--tolerance=" + to_string(aTolerance));

        // This means that ZMiss > 0, which means prediction is activated
        if (aPrediction[0] > 0) {
            argStrings.push_back("--ZMiss=" + to_string(aPrediction[0]));

            // if mspe is activated
            if (aPrediction[1]) {
                argStrings.emplace_back("--mspe");
            }

            // if idw is activated
            if (aPrediction[2]) {
                argStrings.emplace_back("--idw");
            }

            // if fisher is activated
            if (aPrediction[3]) {
                argStrings.emplace_back("--fisher");
            }

            // if mloe_mmom is activated
            if (aPrediction[4]) {
                argStrings.emplace_back("--mloe_mmom");
            }
        }
        // Now, allocate a char** array and copy the arguments into it.
        // Argv ownership is passed to Configurations, and it's responsibility for freeing it.
        char **argv = new char*[argStrings.size()];

        for (size_t i = 0; i < argStrings.size(); ++i) {
            argv[i] = new char[argStrings[i].size() + 1];  // +1 for the null terminator
            std::strcpy(argv[i], argStrings[i].c_str());
        }

        auto configurations = new Configurations();
        configurations->InitializeArguments(argStrings.size(), argv);
        return configurations;
    }

    void R_ExaGeoStatAPI(ExaGeoStatHardware *apHardware, Configurations *apConfigurations) {

        if(apConfigurations->GetPrecision() == exageostat::common::DOUBLE){
            std::unique_ptr<ExaGeoStatData<double>> data;
            exageostat::api::ExaGeoStat<double>::ExaGeoStatLoadData(*apHardware, *apConfigurations, data);
            exageostat::api::ExaGeoStat<double>::ExaGeoStatDataModeling(*apHardware, *apConfigurations, data);
            exageostat::api::ExaGeoStat<double>::ExaGeoStatPrediction(*apHardware, *apConfigurations, data);
        }
        else if(apConfigurations->GetPrecision() == exageostat::common::SINGLE){
            std::unique_ptr<ExaGeoStatData<float>> data;
            exageostat::api::ExaGeoStat<float>::ExaGeoStatLoadData(*apHardware, *apConfigurations, data);
            exageostat::api::ExaGeoStat<float>::ExaGeoStatDataModeling(*apHardware, *apConfigurations, data);
            exageostat::api::ExaGeoStat<float>::ExaGeoStatPrediction(*apHardware, *apConfigurations, data);
        }
    }
}