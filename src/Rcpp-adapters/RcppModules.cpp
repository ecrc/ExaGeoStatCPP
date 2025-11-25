
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file RcppModules.cpp
 * @brief Rcpp module definitions for ExaGeoStatCPP.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author David Helmy
 * @date 2023-01-20
 **/

#include <Rcpp.h>

#include <hardware/ExaGeoStatHardware.hpp>
#include <Rcpp-adapters/FunctionsAdapter.hpp>
#include <data-units/ExaGeoStatData.hpp>

/** Expose C++ class to R to be able to use Wrap and As
 *  Allows C++ to Send and Receive Class object from R
 **/
RCPP_EXPOSED_CLASS(ExaGeoStatHardware)
RCPP_EXPOSED_CLASS_NODECL(ExaGeoStatData<double>)

/** Expose C++ Object With the Given functions **/
RCPP_MODULE(ExaGeoStatCPP) {

    /** ExaGeoStatCPP Class **/
    using namespace Rcpp;

    /** Hardware Class **/
    class_<ExaGeoStatHardware>("Hardware")
        .constructor<std::string, int, int, int, int>()
        .method("finalize_hardware", &ExaGeoStatHardware::FinalizeHardware, 
                "Manually finalize the hardware");

    /** Data Class **/
    class_<ExaGeoStatData<double>>("Data")
        .constructor<int, std::string>();

    function("simulate_data", &exageostat::adapters::R_ExaGeoStatLoadData,
             List::create(
                 _["kernel"], 
                 _["initial_theta"], 
                 _["distance_matrix"] = "euclidean", 
                 _["problem_size"],
                 _["seed"] = static_cast<unsigned int>(time(0)),
                 _["dts"], 
                 _["lts"] = 0,
                 _["dimension"] = "2D", 
                 _["log_path"] = "",
                 _["data_path"] = "", 
                 _["observations_file"] = "", 
                 _["recovery_file"] = ""
             ));

    function("model_data", &exageostat::adapters::R_ExaGeoStatModelData,
             List::create(
                 _["computation"] = "exact", 
                 _["kernel"], 
                 _["distance_matrix"] = "euclidean", 
                 _["lb"], 
                 _["ub"],
                 _["tol"] = 4, 
                 _["mle_itr"], 
                 _["dts"], 
                 _["lts"] = 0, 
                 _["dimension"] = "2D", 
                 _["band"] = 0,
                 _["max_rank"] = 500, 
                 _["acc"] = 0, 
                 _["data"] = R_NilValue, 
                 _["matrix"] = R_NilValue, 
                 _["x"] = R_NilValue,
                 _["y"] = R_NilValue, 
                 _["z"] = R_NilValue
             ));

    function("predict_data", &exageostat::adapters::R_ExaGeoStatPredictData,
             List::create(
                 _["kernel"], 
                 _["distance_matrix"] = "euclidean", 
                 _["estimated_theta"], 
                 _["dts"], 
                 _["lts"] = 0,
                 _["dimension"] = "2D", 
                 _["train_data"], 
                 _["test_data"]
             ));

    function("mloe_mmom", &exageostat::adapters::R_ExaGeoStatMLOE_MMOM,
             List::create(
                 _["kernel"], 
                 _["distance_matrix"] = "euclidean", 
                 _["estimated_theta"], 
                 _["true_theta"],
                 _["dts"], 
                 _["lts"] = 0, 
                 _["dimension"] = "2D", 
                 _["train_data"], 
                 _["test_data"]
             ));

    function("fisher", &exageostat::adapters::R_ExaGeoStatFisher,
             List::create(
                 _["kernel"], 
                 _["distance_matrix"] = "euclidean", 
                 _["estimated_theta"], 
                 _["dts"], 
                 _["lts"] = 0,
                 _["dimension"] = "2D", 
                 _["train_data"], 
                 _["test_data"]
             ));

    function("idw", &exageostat::adapters::R_ExaGeoStatIDW,
             List::create(
                 _["kernel"], 
                 _["distance_matrix"] = "euclidean", 
                 _["estimated_theta"], 
                 _["dts"], 
                 _["lts"] = 0,
                 _["dimension"] = "2D", 
                 _["train_data"], 
                 _["test_data"], 
                 _["test_measurements"]
             ));

    function("mean_trend_removal", &exageostat::adapters::R_ExaGeoStatMeanTrendRemoval,
             List::create(
                 _["lon"],
                 _["startyear"],
                 _["endyear"],
                 _["dts"],
                 _["data_path"],
                 _["forcing_data_path"],
                 _["results_path"],
                 _["starting_theta"],
                 _["lb"],
                 _["ub"],
                 _["tolerance"] = 7,
                 _["max_mle_iterations"] = 30,
                 _["cores"] = 1,
                 _["gpus"] = 0,
                 _["p"] = 1,
                 _["q"] = 1,
                 _["log_path"] = ""
             ));

    function("climate_emulator", &exageostat::adapters::R_ExaGeoStatClimateEmulator,
             List::create(
                 _["N"],
                 _["dts"],
                 _["timeslot"],
                 _["objects_number"],
                 _["data_path"],
                 _["cores"] = 1,
                 _["gpus"] = 0,
                 _["add_diagonal"] = 10.0,
                 _["accuracy"] = 0,
                 _["band_dense_dp"] = 1000,
                 _["hnb"] = 300,
                 _["verbose"] = "detailed",
                 _["log_path"] = ""
             ));

}
