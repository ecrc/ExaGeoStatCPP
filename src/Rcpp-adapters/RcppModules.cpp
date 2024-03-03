
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
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
#include <configurations/Configurations.hpp>
#include <Rcpp-adapters/FunctionsAdapter.hpp>
#include <data-units/ExaGeoStatData.hpp>

/** Expose C++ class to R to be able to use Wrap and As
 *  Allows C++ to Send and Receive Class object from R
 **/
RCPP_EXPOSED_CLASS(ExaGeoStatHardware)
RCPP_EXPOSED_CLASS(Configurations)
RCPP_EXPOSED_CLASS_NODECL(ExaGeoStatData<double>)

/** Expose C++ Object With the Given functions **/
RCPP_MODULE(ExaGeoStatCPP) {

    /** ExaGeoStatCPP Class **/
    using namespace Rcpp;

    /** Hardware Class **/
    class_<ExaGeoStatHardware>("Hardware")
            .constructor<std::string, int, int>()
            .method("finalize_hardware", &ExaGeoStatHardware::FinalizeHardware, "Manually finalize the hardware");

    /** Configurations Class **/
    class_<Configurations>("Configurations")
            .constructor();

    /** Data Class **/
    class_<ExaGeoStatData<double>>("Data")
            .constructor<int, std::string>();

    /** Configurations Function **/
    function("configurations_init", &exageostat::adapters::R_InitializeArguments,
             List::create(_["n"], _["kernel"], _["tile_size"], _["p_q"] = IntegerVector::create(1, 1),
                          _["time_slot"] = 1, _["computation"] = "exact", _["precision"] = "double",
                          _["cores_gpus"] = IntegerVector::create(1, 0), _["band"] = 1, _["max_rank"] = 500,
                          _["iTheta"], _["lb_ub"], _["eTheta"] = IntegerVector::create(-1, -1, -1),
                          _["verbose"] = "standard", _["dimension"] = "2D", _["mle_itr"] = 0, _["tol"] = 4,
                          _["prediction"] = IntegerVector::create(0, 0, 0, 0, 0),
                          _["paths"] = StringVector::create("", "", "", ""), _["save_data"] = 0));

    function("simulate_data", &exageostat::adapters::R_ExaGeoStatLoadData,
             List::create(_["hardware"], _["config"], _["data"]));
    function("model_data", &exageostat::adapters::R_ExaGeoStatModelData,
             List::create(_["hardware"], _["config"], _["data"], _["matrix"] = R_NilValue, _["x"] = R_NilValue,
                          _["y"] = R_NilValue, _["z"] = R_NilValue));
    function("predict_data", &exageostat::adapters::R_ExaGeoStatPredictData,
             List::create(_["hardware"], _["config"], _["data"], _["matrix"] = R_NilValue, _["x"] = R_NilValue,
                          _["y"] = R_NilValue, _["z"] = R_NilValue));

//    function("get_locationsX", &exageostat::adapters::R_GetLocationX, List::create(_["data"]));
//
//    function("get_locationsY", &exageostat::adapters::R_GetLocationY, List::create(_["data"]));
//
    function("get_locationsZ", &exageostat::adapters::R_GetLocationZ, List::create(_["data"]));

//    function("get_Z_measurement_vector", &exageostat::adapters::R_GetDescZValues, List::create(_["data"], _["type"]));

}
