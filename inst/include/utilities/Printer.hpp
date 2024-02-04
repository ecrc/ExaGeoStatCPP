
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file Printer.hpp
 * @version 1.1.0
 * @brief Provides printing functionality for output messages.
 * @details Defines macros for printing messages to the console or R output.
 * @author Mahmoud ElKarargy
 * @author David Helmy
 * @date 2024-01-20
**/

#ifndef EXAGEOSTATCPP_PRINTER_HPP
#define EXAGEOSTATCPP_PRINTER_HPP

#include <sstream>
#ifdef USING_R
//#include <Rcpp.h>
#endif

/**
 * @def EXAGEOSTAT_PRINTER(message)
 * @brief Macro for printing messages.
 * @details Depending on whether the code is compiled with Rcpp support, it prints the message to the console or R output.
 */
#ifdef USING_R
#define EXAGEOSTAT_PRINTER(message) \
    std::cout <<(message);
//    Rcpp::Rcout<<(message);
#else
#define EXAGEOSTAT_PRINTER(message) \
    std::cout<<(message);
#endif

#endif //EXAGEOSTATCPP_PRINTER_HPP
