
# Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file FindPnetCDF.cmake
# @brief A CMake module to locate the Parallel NetCDF library using pkg-config if available.
# @version 2.0.0
# @author Mahmoud ElKarargy
# @date 2024-11-23

# This module defines the following variables:
#  - PnetCDF_FOUND        : True if the library is found
#  - PnetCDF_INCLUDE_DIRS : Path to the include directory
#  - PnetCDF_LIBRARIES    : The library to link against
#  - PnetCDF_VERSION      : Version of the library

find_package(PkgConfig)

if(PkgConfig_FOUND)
    # Use pkg-config to locate PnetCDF
    pkg_check_modules(PnetCDF QUIET pnetcdf)
    if(PnetCDF_FOUND)
        # Assign include directories and libraries from pkg-config
        set(PnetCDF_INCLUDE_DIRS ${PnetCDF_INCLUDE_DIRS} ${PnetCDF_INCLUDEDIR})
        set(PnetCDF_LIBRARIES ${PnetCDF_LIBRARIES})
        set(PnetCDF_VERSION ${PnetCDF_VERSION})
    endif()
endif()

# Fallback if pkg-config failed or did not provide include directories
if(NOT PnetCDF_INCLUDE_DIRS)
    find_path(PnetCDF_INCLUDE_DIR
            NAMES pnetcdf.h
            HINTS ENV PNETCDF_DIR
            PATH_SUFFIXES include
            )
    set(PnetCDF_INCLUDE_DIRS ${PnetCDF_INCLUDE_DIR})
endif()

if(NOT PnetCDF_LIBRARIES)
    find_library(PnetCDF_LIBRARY
            NAMES pnetcdf
            HINTS ENV PNETCDF_DIR
            PATH_SUFFIXES lib
            )
    set(PnetCDF_LIBRARIES ${PnetCDF_LIBRARY})
endif()

if(PnetCDF_INCLUDE_DIRS AND PnetCDF_LIBRARIES)
    set(PnetCDF_FOUND TRUE)
else()
    set(PnetCDF_FOUND FALSE)
endif()

# Detect version from header file if not set
if(PnetCDF_FOUND AND NOT PnetCDF_VERSION)
    file(READ "${PnetCDF_INCLUDE_DIRS}/pnetcdf.h" PNETCDF_HEADER_CONTENTS)
    string(REGEX MATCH "#define PNETCDF_VERSION_MAJOR ([0-9]+)" _major_match "${PNETCDF_HEADER_CONTENTS}")
    string(REGEX MATCH "#define PNETCDF_VERSION_MINOR ([0-9]+)" _minor_match "${PNETCDF_HEADER_CONTENTS}")
    string(REGEX MATCH "#define PNETCDF_VERSION_PATCH ([0-9]+)" _patch_match "${PNETCDF_HEADER_CONTENTS}")
    if(_major_match AND _minor_match AND _patch_match)
        string(REGEX REPLACE ".* ([0-9]+).*" "\\1" PNETCDF_VERSION_MAJOR "${_major_match}")
        string(REGEX REPLACE ".* ([0-9]+).*" "\\1" PNETCDF_VERSION_MINOR "${_minor_match}")
        string(REGEX REPLACE ".* ([0-9]+).*" "\\1" PNETCDF_VERSION_PATCH "${_patch_match}")
        set(PnetCDF_VERSION "${PNETCDF_VERSION_MAJOR}.${PNETCDF_VERSION_MINOR}.${PNETCDF_VERSION_PATCH}")
    endif()
endif()

# Print debug information
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PnetCDF REQUIRED_VARS PnetCDF_INCLUDE_DIRS PnetCDF_LIBRARIES VERSION_VAR PnetCDF_VERSION)

mark_as_advanced(PnetCDF_INCLUDE_DIRS PnetCDF_LIBRARIES)