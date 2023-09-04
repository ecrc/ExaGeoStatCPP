/**
 * @file ChameleonHeaders.hpp
 * @brief This file contains the necessary includes for using the Chameleon library.
 * @version 1.0.0
 * @author Mahmoud ElKarargy
 * @date 2023-08-24
**/

#ifndef EXAGEOSTATCPP_CHAMELEONHEADERS_HPP
#define EXAGEOSTATCPP_CHAMELEONHEADERS_HPP

#ifdef EXAGEOSTAT_USE_CHAMELEON
extern "C" {
#include <chameleon.h>
#include <control/context.h>
#include <control/context.h>
}
#endif

#endif //EXAGEOSTATCPP_CHAMELEONHEADERS_HPP
