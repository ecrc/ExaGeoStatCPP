/**
 * Copyright (c) 2023, King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * MPCR is an R package provided by the STSDS group at KAUST
 *
 **/

#include <Rcpp.h>


using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif


RcppExport SEXP _rcpp_module_boot_ExaGeoStatCPP();

static const R_CallMethodDef CallEntries[] = {
    {"_rcpp_module_boot_ExaGeoStatCPP", (DL_FUNC) &_rcpp_module_boot_ExaGeoStatCPP, 0},
    {NULL, NULL,                                                0}
};

RcppExport void R_init_ExaGeoStatCPP(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
