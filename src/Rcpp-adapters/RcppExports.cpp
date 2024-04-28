
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file RcppExports.cpp
 * @brief Rcpp export definitions for ExaGeoStatCPP module.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author David Helmy
 * @date 2023-01-20
**/

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

/**
 * @brief Rcpp module boot function for ExaGeoStatCPP.
 * @return A SEXP represents the Rcpp module.
 */
RcppExport SEXP _rcpp_module_boot_ExaGeoStatCPP();

/**
 * @brief Array of R function call entries.
 */
static const R_CallMethodDef CallEntries[] = {
        {"_rcpp_module_boot_ExaGeoStatCPP", (DL_FUNC) &_rcpp_module_boot_ExaGeoStatCPP, 0},
        {nullptr,                           nullptr,                                    0}
};

/**
 * @brief R initialization function for ExaGeoStatCPP module.
 * @param dll The DllInfo structure.
 */
RcppExport void R_init_ExaGeoStatCPP(DllInfo *dll) {
    R_registerRoutines(dll, nullptr, CallEntries, nullptr, nullptr);
    R_useDynamicSymbols(dll, FALSE);
}
