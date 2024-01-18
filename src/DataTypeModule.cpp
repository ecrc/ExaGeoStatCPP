/**
 * Copyright (c) 2023, King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * MPCR is an R package provided by the STSDS group at KAUST
 *
 **/
#include <Rcpp.h>
#include <hardware/ExaGeoStatHardware.hpp>

/** Expose C++ class to R to be able to use Wrap and As
 *  Allows C++ to Send and Receive Class object from R
 **/
RCPP_EXPOSED_CLASS(ExaGeoStatHardware)

/** Expose C++ Object With the Given functions **/
RCPP_MODULE(ExaGeoStatCPP) {

    /** MPCR Class **/
    using namespace Rcpp;


    /** Basic Utilities **/
    class_ <exageostat::hardware::ExaGeoStatHardware>("Hardware")
        .constructor <int, int>();

    /** Function that are not masked **/

}
