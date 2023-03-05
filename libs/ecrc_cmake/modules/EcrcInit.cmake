###
#
# @copyright (c) 2009-2014 The University of Tennessee and The University
#                          of Tennessee Research Foundation.
#                          All rights reserved.
# @copyright (c) 2012-2016 Inria. All rights reserved.
# @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
# @copyright (c) 2022 King Abdullah University of Science and Technology (KAUST).
#                     All rights reserved.
#
###
#
#  @file EcrcInit.cmake
#
#  @project ECRC
#  ECRC is a software package provided by:
#     Inria Bordeaux - Sud-Ouest,
#     Univ. of Tennessee,
#     King Abdullah University of Science and Technology
#     Univ. of California Berkeley,
#     Univ. of Colorado Denver.
#
#  @version 0.9.0
#  @author Cedric Castagnede
#  @author Emmanuel Agullo
#  @author Mathieu Faverge
#  @author Florent Pruvost
#  @date 13-07-2012
#
###

# Add the path where we handle our FindFOO.cmake to seek for liraries
if( NOT ECRC_CMAKE_MODULE_PATH )
    message( WARNING "ECRC_CMAKE_MODULE_PATH is not set, assuming default value: \${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules/ecrc/modules" )
    set(ECRC_CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules/ecrc/modules")
endif()
list(APPEND CMAKE_MODULE_PATH ${ECRC_CMAKE_MODULE_PATH}/find)

# This include is required to check symbols of libs in the main CMakeLists.txt
include(CheckFunctionExists)

# This include is required to check defines in headers
include(CheckIncludeFiles)

# To colorize messages
#include(ColorizeMessage)

# To find headers and libs
include(FindHeadersAndLibs)

# Some macros to print status when search for headers and libs
# PrintFindStatus.cmake is in cmake_modules/ecrc/find directory
include(PrintFindStatus)

# Define some auxilary flags
include(AuxilaryFlags)

# Define some variables to et info about ressources
include(Ressources)

option(ECRC_ENABLE_WARNING       "Enable warning messages" OFF)
option(ECRC_ENABLE_COVERAGE      "Enable flags for coverage test" OFF)
#option(ECRC_VERBOSE_FIND_PACKAGE "Add additional messages concerning packages not found" OFF)
#message(STATUS "ECRC_VERBOSE_FIND_PACKAGE is set to OFF, turn it ON to get"
#        "   information about packages not found")

##
## @end file EcrcInit.cmake
##
