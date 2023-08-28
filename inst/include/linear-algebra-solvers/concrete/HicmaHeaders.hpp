/**
 * @file HicmaHeaders.hpp
 * @brief This file contains the necessary includes for using the Chameleon library.
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-08-24
**/

#ifndef EXAGEOSTATCPP_HICMAHEADERS_HPP
#define EXAGEOSTATCPP_HICMAHEADERS_HPP

#ifdef EXAGEOSTAT_USE_HICMA
extern "C" {
#include <hicma_struct.h>
#include <hicma.h>
#include <control/hicma_context.h>
}
#endif
#endif //EXAGEOSTATCPP_HICMAHEADERS_HPP
