
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file TestPrediction.cpp
 * @brief Unit tests for the TestPrediction class in the ExaGeoStat software package.
 * @details This file contains Catch2 unit tests that validate the functionality of the TestPrediction class
 * @version 1.0.0
 * @author Mahmoud ElKarargy
 * @date 2024-1-18
**/

#include <catch2/catch_all.hpp>

#include <prediction/Prediction.hpp>

using namespace std;

using namespace exageostat::common;
using namespace exageostat::configurations;
using namespace exageostat::hardware;
using namespace exageostat::prediction;

void TEST_PREDICTION_MISSING_DATA() {

    //Init configuration
    Configurations configurations;
    configurations.SetUnknownObservationsNb(4);
    int N = 16;
    configurations.SetProblemSize(N);
    configurations.SetKernelName("UnivariateMaternStationary");
    int dts = 8;

    configurations.SetDenseTileSize(dts);
    configurations.SetComputation(EXACT_DENSE);
    configurations.SetMaxMleIterations(3);
    configurations.SetTolerance(pow(10, -4));

    vector<double> lb{0.1, 0.1, 0.1};
    configurations.SetLowerBounds(lb);
    configurations.SetStartingTheta(lb);
    vector<double> ub{5, 5, 5};
    configurations.SetUpperBounds(ub);
    vector<double> initial_theta{1, 0.1, 0.5};
    configurations.SetInitialTheta(initial_theta);


    auto hardware = ExaGeoStatHardware(EXACT_DENSE, configurations.GetCoresNumber(), configurations.GetGPUsNumbers());
    std::unique_ptr<exageostat::dataunits::ExaGeoStatData<double>> data = std::make_unique<exageostat::dataunits::ExaGeoStatData<double>>(
            configurations.GetProblemSize(), configurations.GetDimension());

    auto *z_matrix = new double[N]{-1.272336140360187606, -2.590699695867695773, 0.512142584178685967,
                                   -0.163880452049749520, 0.313503633252489700, -1.474410682226017677,
                                   0.161705025505231914, 0.623389205185149065, -1.341858445399783495,
                                   -1.054282062428600009, -1.669383221392507943, 0.219170645803740793,
                                   0.971213790000161170, 0.538973474182433021, -0.752828466476077041,
                                   0.290822066007430102};
    //creating locations x and y.
    auto *location_x = new double[N]{0.193041886015106440, 0.330556191348134576, 0.181612878614480805,
                                     0.370473792629892440, 0.652140077821011688, 0.806332494087129037,
                                     0.553322652018005678, 0.800961318379491916, 0.207324330510414295,
                                     0.347951476310368490, 0.092042420080872822, 0.465445944914930965,
                                     0.528267338063630132, 0.974792095826657490, 0.552452887769893985,
                                     0.877592126344701295};
    auto *location_y = new double[N]{0.103883421072709245, 0.135790035858701447, 0.434683756771190977,
                                     0.400778210116731537, 0.168459601739528508, 0.105195696955825133,
                                     0.396398870832379624, 0.296757457846952011, 0.564507515068284116,
                                     0.627679865720607300, 0.928648813611047563, 0.958236057068741931,
                                     0.573571374074921758, 0.568657969024185528, 0.935835812924391552,
                                     0.942824444953078489};

    data->GetLocations()->SetLocationX(*location_x, N);
    data->GetLocations()->SetLocationY(*location_y, N);

    SECTION("Test Prediction - MSPE ")
    {
        configurations.SetIsIDW(false);
        configurations.SetIsMLOEMMOM(false);
        configurations.SetIsFisher(false);
        configurations.SetIsMSPE(true);

        vector<double> estimated_theta{0.9, 0.09, 0.4};
        configurations.SetEstimatedTheta(estimated_theta);

        Prediction<double> predictor;

        // Register and create a kernel object
        exageostat::kernels::Kernel<double> *pKernel = exageostat::plugins::PluginRegistry<exageostat::kernels::Kernel<double>>::Create(
                configurations.GetKernelName(),
                configurations.GetTimeSlot());

        // Add the data prediction arguments.
        configurations.InitializeDataPredictionArguments();
        predictor.PredictMissingData(hardware, data, configurations, z_matrix, *pKernel);

        REQUIRE(exageostat::results::Results::GetInstance()->GetMSPEError()== Catch::Approx(0.552448));
        delete pKernel;
    }
    SECTION("Test Prediction - IDW ") {
        configurations.SetIsMLOEMMOM(false);
        configurations.SetIsFisher(false);
        configurations.SetIsMSPE(false);
        configurations.SetIsIDW(true);

        vector<double> estimated_theta{0.9, 0.09, 0.4};
        configurations.SetEstimatedTheta(estimated_theta);

        Prediction<double> predictor;

        std::vector<double> idw_error={ 1.18856255, 1.25725881, 1.11986628 };

        // Register and create a kernel object
        exageostat::kernels::Kernel<double> *pKernel = exageostat::plugins::PluginRegistry<exageostat::kernels::Kernel<double>>::Create(
                configurations.GetKernelName(),
                configurations.GetTimeSlot());
        // Add the data prediction arguments.
        configurations.InitializeDataPredictionArguments();
        predictor.PredictMissingData(hardware, data, configurations, z_matrix, *pKernel);
        for(int i =0;i<3;i++){
            REQUIRE(exageostat::results::Results::GetInstance()->GetIDWError()[i] == Catch::Approx(idw_error[i]));
        }
        delete pKernel;
    }
    SECTION("Test Prediction - MLOE_MMOM ") {
        configurations.SetIsMSPE(false);
        configurations.SetIsIDW(false);
        configurations.SetIsMLOEMMOM(true);
        configurations.SetIsFisher(false);

        vector<double> estimated_theta{0.9, 0.09, 0.4};
        configurations.SetEstimatedTheta(estimated_theta);

        Prediction<double> predictor;
        // Register and create a kernel object
        exageostat::kernels::Kernel<double> *pKernel = exageostat::plugins::PluginRegistry<exageostat::kernels::Kernel<double>>::Create(
                configurations.GetKernelName(),
                configurations.GetTimeSlot());
        // Add the data prediction arguments.
        configurations.InitializeDataPredictionArguments();
        predictor.PredictMissingData(hardware, data, configurations, z_matrix, *pKernel);
        REQUIRE(exageostat::results::Results::GetInstance()->GetMLOE() == Catch::Approx(0.004467).margin(0.001));
        REQUIRE(exageostat::results::Results::GetInstance()->GetMMOM() == Catch::Approx(-0.0812376).margin(0.001));


        vector<double> new_estimated_theta1{1, 0.1, 0.5};
        configurations.SetEstimatedTheta(new_estimated_theta1);
        // Add the data prediction arguments.
        configurations.InitializeDataPredictionArguments();
        predictor.PredictMissingData(hardware, data, configurations, z_matrix, *pKernel);
        REQUIRE(exageostat::results::Results::GetInstance()->GetMLOE() == Catch::Approx(0).margin(0.001));
        REQUIRE(exageostat::results::Results::GetInstance()->GetMMOM() == Catch::Approx(0).margin(0.001));
        delete pKernel;
    }

    SECTION("Test Prediction - FISHER") {
        configurations.SetIsMSPE(false);
        configurations.SetIsIDW(false);
        configurations.SetIsMLOEMMOM(false);
        configurations.SetIsFisher(true);

        vector<double> new_estimated_theta{0.9, 0.09, 0.4};
        configurations.SetEstimatedTheta(new_estimated_theta);
        Prediction<double> predictor;
        // Register and create a kernel object
        exageostat::kernels::Kernel<double> *pKernel = exageostat::plugins::PluginRegistry<exageostat::kernels::Kernel<double>>::Create(
                configurations.GetKernelName(),
                configurations.GetTimeSlot());
        // Add the data prediction arguments.
        configurations.InitializeDataPredictionArguments();
        predictor.PredictMissingData(hardware, data, configurations, z_matrix, *pKernel);

        REQUIRE(exageostat::results::Results::GetInstance()->GetFisher00() == Catch::Approx(0.104589));
        REQUIRE(exageostat::results::Results::GetInstance()->GetFisher11() == Catch::Approx(0.187355));
        REQUIRE(exageostat::results::Results::GetInstance()->GetFisher22() == Catch::Approx(10.556483));
        delete pKernel;
    }
    SECTION("Test Prediction - Exception"){
        vector<double> new_estimated_theta{-1, -1, -1};
        configurations.SetEstimatedTheta(new_estimated_theta);
        Prediction<double> predictor;
        // Register and create a kernel object
        configurations.SetIsMLOEMMOM(true);
        exageostat::kernels::Kernel<double> *pKernel = exageostat::plugins::PluginRegistry<exageostat::kernels::Kernel<double>>::Create(
                configurations.GetKernelName(),
                configurations.GetTimeSlot());
        REQUIRE_THROWS(predictor.PredictMissingData(hardware, data, configurations, z_matrix, *pKernel));
        delete pKernel;
    }
    delete[] location_x;
    delete[] location_y;
    delete[] z_matrix;
}

TEST_CASE("Test Predictions") {
    TEST_PREDICTION_MISSING_DATA();
}