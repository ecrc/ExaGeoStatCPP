
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

/** Expose C++ class to R to be able to use Wrap and As
 *  Allows C++ to Send and Receive Class object from R
 **/

RCPP_EXPOSED_CLASS(ExaGeoStatHardware)
RCPP_EXPOSED_CLASS(Configurations)

/** Expose C++ Object With the Given functions **/
RCPP_MODULE(ExaGeoStatCPP) {

    /** ExaGeoStatCPP Class **/
    using namespace Rcpp;

    /** Hardware Class **/
    class_<ExaGeoStatHardware>("Hardware")
            .constructor<std::string, int, int>();

    /** Configurations Class **/
    class_<Configurations>("Configurations")
            .constructor();

    /** Configurations Function **/
    function("configurations_init", &exageostat::adapters::R_InitializeArguments,
             List::create(_["n"], _["kernel"], _["tile_size"], _["p_q"] = IntegerVector::create(1, 1), _["time_slot"] = 1,
                          _["computation"] = "exact", _["precision"] = "double", _["cores_gpus"] = IntegerVector::create(1, 0),
                          _["band"] = 0, _["max_rank"] = 1, _["iTheta"], _["lb_ub"],
                          _["eTheta"]=IntegerVector::create(-1,-1,-1), _["verbose"] = "standard", _["dimension"] = "2D",
                          _["mle_itr"], _["tol"] = 4, _["prediction"] = IntegerVector::create(0, 0, 0, 0, 0)));

    function("exageostat", &exageostat::adapters::R_ExaGeoStatAPI, List::create(_["hardware"], _["config"]));
}
