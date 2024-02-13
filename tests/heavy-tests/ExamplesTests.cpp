
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file TestExaGeoStatApi.cpp
 * @brief Test suite for the ExaGeoStat APIs data generation functionality.
 * @version 1.0.1
 * @author Mahmoud ElKarargy
 * @date 2023-12-26
**/

#include <iostream>
#include <catch2/catch_all.hpp>

#include <api/ExaGeoStat.hpp>

using namespace std;
using namespace std::filesystem;

using namespace exageostat::configurations;
using namespace exageostat::dataunits;
using namespace exageostat::api;

TEST_CASE("EXAMPLES") {

    // Specify the examples path
    path currentPath = current_path().parent_path().parent_path() / "examples";

    SECTION("Configuration"){
        path configurations_example = currentPath / "configurations/Example_Configurations_Setup ";
        string arguments_string = "--N=10 --dts=8 --kernel=univariate_matern_stationary";
        cout << "Running Configurations example with arguments: " + arguments_string << endl << flush;

        std::string fullCommand = std::string(configurations_example) + " " + arguments_string;
        if (std::system(fullCommand.c_str())) {
            throw runtime_error("This test failed " + fullCommand);
        }
    }
    SECTION("Synthetic Data Generation"){
        path synthetic_data_example = currentPath / "data-generators/Example_Synthetic_Data_Generation";
        string arguments_string = "--N=16 --dts=8 --kernel=univariate_matern_stationary --dimension=2d --computation=exact --itheta=1:0.1:0.5 --ub=5:5:5 --lb=0.1:0.1:0.1";
        cout << "Running synthetic locations example with arguments: " + arguments_string << endl << flush;

        std::string fullCommand = std::string(synthetic_data_example) + " " + arguments_string;
        if (std::system(fullCommand.c_str())) {
            throw runtime_error("This test failed " + fullCommand);
        }
    }
    currentPath = currentPath / "end-to-end";
    SECTION("Data Generation Module"){
        path data_generation_example = currentPath / "Example_Data_Generation";
        string arguments_string = "--N=32 --dts=8 --kernel=univariate_matern_stationary --dimension=3D --computation=diagonal_approx --itheta=1:0.1:0.5 --ub=5:5:5 --lb=0.1:0.1:0.1";
        cout << "Running Data Generation Module example with arguments: " + arguments_string << endl << flush;

        std::string fullCommand = std::string(data_generation_example) + " " + arguments_string;
        if (std::system(fullCommand.c_str())) {
            throw runtime_error("This test failed " + fullCommand);
        }
    }
    SECTION("Data Modeling Module"){
        path data_modeling_example = currentPath / "Example_Data_Modeling";
#ifdef USE_HICMA
        string arguments_string = "--N=16 --dts=8 --lts=8 --kernel=univariate_matern_stationary --dimension=2D --computation=tlr --itheta=1:0.1:0.5 --ub=5:5:5 --lb=0.1:0.1:0.1 --max_mle_iterations=3 --tolerance=10 --max_rank=500";
#else
        string arguments_string = "--N=16 --dts=8 --kernel=univariate_matern_stationary --computation=exact --itheta=1:0.1:0.5 --ub=5:5:5 --lb=0.1:0.1:0.1 --max_mle_iterations=3 --tolerance=10";
#endif
        cout << "Running Data Modeling example with arguments: " + arguments_string << endl << flush;
        std::string fullCommand = std::string(data_modeling_example) + " " + arguments_string;
        if (std::system(fullCommand.c_str())) {
            throw runtime_error("This test failed " + fullCommand);
        }
    }
    SECTION("Data Prediction Module"){
        path data_prediction_example = currentPath / "Example_Data_Prediction";

        string arguments_string = "--N=16 --dts=8 --kernel=univariate_matern_stationary --computation=exact --itheta=1:0.1:0.5 --ub=5:5:5 --lb=0.1:0.1:0.1 --Zmiss=10 --mspe --fisher --mloe_mmom --etheta=1.1:0.1:0.5";
        cout << "Running Data Prediction example with arguments: " + arguments_string << endl << flush;

        std::string fullCommand = std::string(data_prediction_example) + " " + arguments_string;
        if (std::system(fullCommand.c_str())) {
            throw runtime_error("This test failed " + fullCommand);
        }
    }
    SECTION("Data Generation and Modeling"){
        path data_generation_modeling_example = currentPath / "Example_Data_Generation_and_Modeling";

        string arguments_string = "--N=16 --dts=8 --kernel=univariate_matern_stationary --computation=exact --itheta=1:0.1:0.5 --ub=5:5:5 --lb=0.1:0.1:0.1 --max_mle_iterations=2 --tolerance=4 --etheta=1.1:?:0.5";
        cout << "Running Data data generation and modeling example with arguments: " + arguments_string << endl << flush;

        std::string fullCommand = std::string(data_generation_modeling_example) + " " + arguments_string;
        if (std::system(fullCommand.c_str())) {
            throw runtime_error("This test failed " + fullCommand);
        }
    }
}

