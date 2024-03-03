/**
 * @file HicmaHeaders.hpp
 * @brief This file contains the necessary includes for using the Chameleon library.
 * @version 1.0.1
 * @author Mahmoud ElKarargy
 * @date 2023-08-24
**/

#ifndef EXAGEOSTATCPP_HICMAHEADERS_HPP
#define EXAGEOSTATCPP_HICMAHEADERS_HPP

#ifdef USE_HICMA
extern "C" {
#include <hicma_struct.h>
#include <starsh-spatial.h>
#include <hicma_z.h>
#include "hicma/misc/auxdescutil.h"
#include <hicma_d.h>
#include <control/hicma_descriptor.h>
#include <control/hicma_context.h>
#include "hicma/misc/auxcompute_z.h"
}
#endif
#endif //EXAGEOSTATCPP_HICMAHEADERS_HPP
