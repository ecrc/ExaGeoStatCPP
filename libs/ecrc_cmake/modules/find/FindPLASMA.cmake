###
#
# @copyright (c) 2009-2014 The University of Tennessee and The University
#                          of Tennessee Research Foundation.
#                          All rights reserved.
# @copyright (c) 2012-2014 Inria. All rights reserved.
# @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
# @copyright (c) 2022 King Abdullah University of Science and Technology (KAUST).
#                     All rights reserved.
#
#
###
#
# - Find PLASMA include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(PLASMA
#               [REQUIRED]             # Fail with error if plasma is not found
#               [COMPONENTS <comp1> <comp2> ...] # dependencies
#              )
#
#  PLASMA depends on the following libraries:
#   - LAPACK
#   - LAPACKE
#   - BLAS
#   - CBLAS
#
#  COMPONENTS are optional libraries PLASMA could be linked with,
#  Use it to drive detection of a specific compilation chain
#  COMPONENTS can be some of the following:
#   - no components are available for now: maybe PLASMA in the future?
#
# Results are reported in variables:
#  PLASMA_FOUND            - True if headers and requested libraries were found
#  PLASMA_LINKER_FLAGS     - list of required linker flags (excluding -l and -L)
#  PLASMA_INCLUDE_DIRS     - plasma include directories
#  PLASMA_LIBRARY_DIRS     - Link directories for plasma libraries
#  PLASMA_LIBRARIES        - plasma libraries
#  PLASMA_INCLUDE_DIRS_DEP - plasma + dependencies include directories
#  PLASMA_LIBRARY_DIRS_DEP - plasma + dependencies link directories
#  PLASMA_LIBRARIES_DEP    - plasma libraries + dependencies
#
# The user can give specific paths where to find the libraries adding cmake
# options at configure (ex: cmake path/to/project -DPLASMA_DIR=path/to/plasma):
#  PLASMA_DIR              - Where to find the base directory of plasma
#  PLASMA_INCDIR           - Where to find the header files
#  PLASMA_LIBDIR           - Where to find the library files
# The module can also look for the following environment variables if paths
# are not given as cmake variable: PLASMA_DIR, PLASMA_INCDIR, PLASMA_LIBDIR
#
#=============================================================================
# Copyright 2012-2013 Inria
# Copyright 2012-2013 Emmanuel Agullo
# Copyright 2012-2013 Mathieu Faverge
# Copyright 2012      Cedric Castagnede
# Copyright 2013      Florent Pruvost
# Copyright 2017      Eduardo Gonzalez
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file MORSE-Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of Morse, substitute the full
#  License text for the above reference.)


if(NOT PLASMA_FOUND)
    set(PLASMA_DIR "" CACHE PATH "Installation directory of PLASMA library")
    if (NOT PLASMA_FIND_QUIETLY)
        message(STATUS "A cache variable, namely PLASMA_DIR, has been set to specify the install directory of PLASMA")
    endif()
endif(NOT PLASMA_FOUND)

# PLASMA depends on LAPACKE anyway, try to find it
if (NOT LAPACKE_FOUND)
    if(PLASMA_FIND_REQUIRED)
        find_package(LAPACKE REQUIRED)
    else()
        find_package(LAPACKE)
    endif()
endif()
# PLASMA depends on CBLAS anyway, try to find it
if (NOT CBLAS_FOUND)
    if(PLASMA_FIND_REQUIRED)
        find_package(CBLAS REQUIRED)
    else()
        find_package(CBLAS)
    endif()
endif()
# BLAS and LAPACK are searched by CBLAS and LAPACKE


set(ENV_PLASMA_DIR "$ENV{PLASMA_DIR}")
set(ENV_PLASMA_INCDIR "$ENV{PLASMA_INCDIR}")
set(ENV_PLASMA_LIBDIR "$ENV{PLASMA_LIBDIR}")
set(PLASMA_GIVEN_BY_USER "FALSE")
if ( PLASMA_DIR OR ( PLASMA_INCDIR AND PLASMA_LIBDIR) OR ENV_PLASMA_DIR OR (ENV_PLASMA_INCDIR AND ENV_PLASMA_LIBDIR) )
    set(PLASMA_GIVEN_BY_USER "TRUE")
endif()

# Optionally use pkg-config to detect include/library dirs (if pkg-config is available)
# -------------------------------------------------------------------------------------
include(FindPkgConfig)
find_package(PkgConfig QUIET)
if(PKG_CONFIG_EXECUTABLE AND NOT PLASMA_GIVEN_BY_USER)

    pkg_search_module(PLASMA plasma)
    if (NOT PLASMA_FIND_QUIETLY)
        if (PLASMA_FOUND AND PLASMA_LIBRARIES)
            message(STATUS "Looking for PLASMA - found using PkgConfig")
            #if(NOT PLASMA_INCLUDE_DIRS)
            #    message("${Magenta}PLASMA_INCLUDE_DIRS is empty using PkgConfig."
            #        "Perhaps the path to plasma headers is already present in your"
            #        "C(PLUS)_INCLUDE_PATH environment variable.${ColourReset}")
            #endif()
        else()
            message("${Magenta}Looking for PLASMA - not found using PkgConfig. "
                "Perhaps you should add the directory containing plasma.pc "
                "to the PKG_CONFIG_PATH environment variable.${ColourReset}")
        endif()
    endif()

    if (PLASMA_FIND_VERSION_EXACT)
        if( NOT (PLASMA_FIND_VERSION_MAJOR STREQUAL PLASMA_VERSION_MAJOR) OR
            NOT (PLASMA_FIND_VERSION_MINOR STREQUAL PLASMA_VERSION_MINOR) )
            if(NOT PLASMA_FIND_QUIETLY)
                message(FATAL_ERROR
                        "PLASMA version found is ${PLASMA_VERSION_STRING} "
                        "when required is ${PLASMA_FIND_VERSION}")
            endif()
        endif()
    else()
        # if the version found is older than the required then error
        if( (PLASMA_FIND_VERSION_MAJOR STRGREATER PLASMA_VERSION_MAJOR) OR
            (PLASMA_FIND_VERSION_MINOR STRGREATER PLASMA_VERSION_MINOR) )
            if(NOT PLASMA_FIND_QUIETLY)
                message(FATAL_ERROR
                        "PLASMA version found is ${PLASMA_VERSION_STRING} "
                        "when required is ${PLASMA_FIND_VERSION} or newer")
            endif()
        endif()
    endif()

    # if pkg-config is used: these variables are empty
    # the pkg_search_module call will set the following:
    # PLASMA_LDFLAGS: all required linker flags
    # PLASMA_CFLAGS:  all required cflags
    set(PLASMA_INCLUDE_DIRS_DEP "")
    set(PLASMA_LIBRARY_DIRS_DEP "")
    set(PLASMA_LIBRARIES_DEP "")
    # replace it anyway: we should update it with dependencies given by pkg-config
    set(PLASMA_INCLUDE_DIRS_DEP "${PLASMA_INCLUDE_DIRS}")
    set(PLASMA_LIBRARY_DIRS_DEP "${PLASMA_LIBRARY_DIRS}")
    set(PLASMA_LIBRARIES_DEP "${PLASMA_LIBRARIES}")

endif(PKG_CONFIG_EXECUTABLE AND NOT PLASMA_GIVEN_BY_USER)

# if PLASMA is not found using pkg-config
if( (NOT PKG_CONFIG_EXECUTABLE) OR (PKG_CONFIG_EXECUTABLE AND NOT PLASMA_FOUND) OR (PLASMA_GIVEN_BY_USER) )

    if (NOT PLASMA_FIND_QUIETLY)
        message(STATUS "Looking for PLASMA - PkgConfig not used")
    endif()

    # Looking for include
    # -------------------

    # Add system include paths to search include
    # ------------------------------------------
    unset(_inc_env)
    set(ENV_PLASMA_DIR "$ENV{PLASMA_DIR}")
    set(ENV_PLASMA_INCDIR "$ENV{PLASMA_INCDIR}")
    if(ENV_PLASMA_INCDIR)
        list(APPEND _inc_env "${ENV_PLASMA_INCDIR}")
    elseif(ENV_PLASMA_DIR)
        list(APPEND _inc_env "${ENV_PLASMA_DIR}")
        list(APPEND _inc_env "${ENV_PLASMA_DIR}/include")
        list(APPEND _inc_env "${ENV_PLASMA_DIR}/include/plasma")
    else()
        if(WIN32)
            string(REPLACE ":" ";" _inc_env "$ENV{INCLUDE}")
        else()
            string(REPLACE ":" ";" _path_env "$ENV{INCLUDE}")
            list(APPEND _inc_env "${_path_env}")
            string(REPLACE ":" ";" _path_env "$ENV{C_INCLUDE_PATH}")
            list(APPEND _inc_env "${_path_env}")
            string(REPLACE ":" ";" _path_env "$ENV{CPATH}")
            list(APPEND _inc_env "${_path_env}")
            string(REPLACE ":" ";" _path_env "$ENV{INCLUDE_PATH}")
            list(APPEND _inc_env "${_path_env}")
        endif()
    endif()
    list(APPEND _inc_env "${CMAKE_PLATFORM_IMPLICIT_INCLUDE_DIRECTORIES}")
    list(APPEND _inc_env "${CMAKE_C_IMPLICIT_INCLUDE_DIRECTORIES}")
    list(REMOVE_DUPLICATES _inc_env)


    # Try to find the plasma header in the given paths
    # -------------------------------------------------
    # call cmake macro to find the header path
    if(PLASMA_INCDIR)
        set(PLASMA_plasma.h_DIRS "PLASMA_plasma.h_DIRS-NOTFOUND")
        find_path(PLASMA_plasma.h_DIRS
          NAMES plasma.h
          HINTS ${PLASMA_INCDIR})
    else()
        if(PLASMA_DIR)
            set(PLASMA_plasma.h_DIRS "PLASMA_plasma.h_DIRS-NOTFOUND")
            find_path(PLASMA_plasma.h_DIRS
              NAMES plasma.h
              HINTS ${PLASMA_DIR}
              PATH_SUFFIXES "include" "include/plasma")
        else()
            set(PLASMA_plasma.h_DIRS "PLASMA_plasma.h_DIRS-NOTFOUND")
            find_path(PLASMA_plasma.h_DIRS
              NAMES plasma.h
              HINTS ${_inc_env})
        endif()
    endif()
    mark_as_advanced(PLASMA_plasma.h_DIRS)

    # If found, add path to cmake variable
    # ------------------------------------
    if (PLASMA_plasma.h_DIRS)
        set(PLASMA_INCLUDE_DIRS "${PLASMA_plasma.h_DIRS}")
    else ()
        set(PLASMA_INCLUDE_DIRS "PLASMA_INCLUDE_DIRS-NOTFOUND")
        if(NOT PLASMA_FIND_QUIETLY)
            message(STATUS "Looking for plasma -- plasma.h not found")
        endif()
    endif()


    # Looking for lib
    # ---------------

    # Add system library paths to search lib
    # --------------------------------------
    unset(_lib_env)
    set(ENV_PLASMA_LIBDIR "$ENV{PLASMA_LIBDIR}")
    if(ENV_PLASMA_LIBDIR)
        list(APPEND _lib_env "${ENV_PLASMA_LIBDIR}")
    elseif(ENV_PLASMA_DIR)
        list(APPEND _lib_env "${ENV_PLASMA_DIR}")
        list(APPEND _lib_env "${ENV_PLASMA_DIR}/lib")
    else()
        if(WIN32)
            string(REPLACE ":" ";" _lib_env "$ENV{LIB}")
        else()
            if(APPLE)
                string(REPLACE ":" ";" _lib_env "$ENV{DYLD_LIBRARY_PATH}")
            else()
                string(REPLACE ":" ";" _lib_env "$ENV{LD_LIBRARY_PATH}")
            endif()
            list(APPEND _lib_env "${CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES}")
            list(APPEND _lib_env "${CMAKE_C_IMPLICIT_LINK_DIRECTORIES}")
        endif()
    endif()
    list(REMOVE_DUPLICATES _lib_env)

    # Try to find the plasma lib in the given paths
    # ----------------------------------------------

    # call cmake macro to find the lib path
    if(PLASMA_LIBDIR)
        set(PLASMA_plasma_LIBRARY "PLASMA_plasma_LIBRARY-NOTFOUND")
        find_library(PLASMA_plasma_LIBRARY
            NAMES plasma
            HINTS ${PLASMA_LIBDIR})
    else()
        if(PLASMA_DIR)
            set(PLASMA_plasma_LIBRARY "PLASMA_plasma_LIBRARY-NOTFOUND")
            find_library(PLASMA_plasma_LIBRARY
                NAMES plasma
                HINTS ${PLASMA_DIR}
                PATH_SUFFIXES lib lib32 lib64)
        else()
            set(PLASMA_plasma_LIBRARY "PLASMA_plasma_LIBRARY-NOTFOUND")
            find_library(PLASMA_plasma_LIBRARY
                NAMES plasma
                HINTS ${_lib_env})
        endif()
    endif()
    mark_as_advanced(PLASMA_plasma_LIBRARY)

    # If found, add path to cmake variable
    # ------------------------------------
    if (PLASMA_plasma_LIBRARY)
        get_filename_component(plasma_lib_path "${PLASMA_plasma_LIBRARY}" PATH)
        # set cmake variables
        set(PLASMA_LIBRARIES    "${PLASMA_plasma_LIBRARY}")
        set(PLASMA_LIBRARY_DIRS "${plasma_lib_path}")
    else ()
        set(PLASMA_LIBRARIES    "PLASMA_LIBRARIES-NOTFOUND")
        set(PLASMA_LIBRARY_DIRS "PLASMA_LIBRARY_DIRS-NOTFOUND")
        if(NOT PLASMA_FIND_QUIETLY)
            message(STATUS "Looking for plasma -- lib plasma not found")
        endif()
    endif ()

    # check a function to validate the find
    if (PLASMA_LIBRARIES)

        set(REQUIRED_LDFLAGS)
        set(REQUIRED_INCDIRS)
        set(REQUIRED_LIBDIRS)
        set(REQUIRED_LIBS)

        # PLASMA
        if (PLASMA_INCLUDE_DIRS)
            set(REQUIRED_INCDIRS "${PLASMA_INCLUDE_DIRS}")
        endif()
        if (PLASMA_LIBRARY_DIRS)
            set(REQUIRED_LIBDIRS "${PLASMA_LIBRARY_DIRS}")
        endif()
        set(REQUIRED_LIBS "${PLASMA_LIBRARIES}")
        # CBLAS
        if (CBLAS_INCLUDE_DIRS_DEP)
            list(APPEND REQUIRED_INCDIRS "${CBLAS_INCLUDE_DIRS_DEP}")
        elseif (CBLAS_INCLUDE_DIRS)
            list(APPEND REQUIRED_INCDIRS "${CBLAS_INCLUDE_DIRS}")
        endif()
        if(CBLAS_LIBRARY_DIRS_DEP)
            list(APPEND REQUIRED_LIBDIRS "${CBLAS_LIBRARY_DIRS_DEP}")
        elseif(CBLAS_LIBRARY_DIRS)
            list(APPEND REQUIRED_LIBDIRS "${CBLAS_LIBRARY_DIRS}")
        endif()
        if (CBLAS_LIBRARIES_DEP)
            list(APPEND REQUIRED_LIBS "${CBLAS_LIBRARIES_DEP}")
        elseif(CBLAS_LIBRARIES)
            list(APPEND REQUIRED_LIBS "${CBLAS_LIBRARIES}")
        endif()
        if (BLAS_LINKER_FLAGS)
            list(APPEND REQUIRED_LDFLAGS "${BLAS_LINKER_FLAGS}")
        endif()
        # LAPACK
        if (LAPACK_INCLUDE_DIRS)
            list(APPEND REQUIRED_INCDIRS "${LAPACK_INCLUDE_DIRS}")
        endif()
        if(LAPACK_LIBRARY_DIRS)
            list(APPEND REQUIRED_LIBDIRS "${LAPACK_LIBRARY_DIRS}")
        endif()
        list(APPEND REQUIRED_LIBS "${LAPACK_LIBRARIES}")
        if (LAPACK_LINKER_FLAGS)
            list(APPEND REQUIRED_LDFLAGS "${LAPACK_LINKER_FLAGS}")
        endif()

        # set required libraries for link
        set(CMAKE_REQUIRED_INCLUDES "${REQUIRED_INCDIRS}")
        set(CMAKE_REQUIRED_LIBRARIES)
        list(APPEND CMAKE_REQUIRED_LIBRARIES "${REQUIRED_LDFLAGS}")
        foreach(lib_dir ${REQUIRED_LIBDIRS})
            list(APPEND CMAKE_REQUIRED_LIBRARIES "-L${lib_dir}")
        endforeach()
        list(APPEND CMAKE_REQUIRED_LIBRARIES "${REQUIRED_LIBS}")
        string(REGEX REPLACE "^ -" "-" CMAKE_REQUIRED_LIBRARIES "${CMAKE_REQUIRED_LIBRARIES}")

        # test link
        unset(PLASMA_WORKS CACHE)
        include(CheckFunctionExists)
        check_function_exists(plasma_dgetrf PLASMA_WORKS)
        mark_as_advanced(PLASMA_WORKS)

        if(PLASMA_WORKS)
            # save link with dependencies
            set(PLASMA_LIBRARIES_DEP    "${REQUIRED_LIBS}")
            set(PLASMA_LIBRARY_DIRS_DEP "${REQUIRED_LIBDIRS}")
            set(PLASMA_INCLUDE_DIRS_DEP "${REQUIRED_INCDIRS}")
            set(PLASMA_LINKER_FLAGS     "${REQUIRED_LDFLAGS}")
            list(REMOVE_DUPLICATES PLASMA_LIBRARY_DIRS_DEP)
            list(REMOVE_DUPLICATES PLASMA_INCLUDE_DIRS_DEP)
            list(REMOVE_DUPLICATES PLASMA_LINKER_FLAGS)
        else()
            if(NOT PLASMA_FIND_QUIETLY)
                message(STATUS "Looking for plasma : test of plasma_dgetrf with
                plasma, cblas, and lapack libraries fails")
                message(STATUS "CMAKE_REQUIRED_LIBRARIES: ${CMAKE_REQUIRED_LIBRARIES}")
                message(STATUS "CMAKE_REQUIRED_INCLUDES: ${CMAKE_REQUIRED_INCLUDES}")
                message(STATUS "Check in CMakeFiles/CMakeError.log to figure out why it fails")
            endif()
        endif()
        set(CMAKE_REQUIRED_INCLUDES)
        set(CMAKE_REQUIRED_FLAGS)
        set(CMAKE_REQUIRED_LIBRARIES)
    endif(PLASMA_LIBRARIES)

endif( (NOT PKG_CONFIG_EXECUTABLE) OR (PKG_CONFIG_EXECUTABLE AND NOT PLASMA_FOUND) OR (PLASMA_GIVEN_BY_USER) )

if (PLASMA_LIBRARIES)
    if (PLASMA_LIBRARY_DIRS)
        foreach(dir ${PLASMA_LIBRARY_DIRS})
            if ("${dir}" MATCHES "plasma")
                set(first_lib_path "${dir}")
            endif()
        endforeach()
    else()
        list(GET PLASMA_LIBRARIES 0 first_lib)
        get_filename_component(first_lib_path "${first_lib}" PATH)
    endif()
    if (${first_lib_path} MATCHES "/lib(32|64)?$")
        string(REGEX REPLACE "/lib(32|64)?$" "" not_cached_dir "${first_lib_path}")
        set(PLASMA_DIR_FOUND "${not_cached_dir}" CACHE PATH "Installation directory of PLASMA library" FORCE)
    else()
        set(PLASMA_DIR_FOUND "${first_lib_path}" CACHE PATH "Installation directory of PLASMA library" FORCE)
    endif()
endif()

# check that PLASMA has been found
# -------------------------------
include(FindPackageHandleStandardArgs)
if (PKG_CONFIG_EXECUTABLE AND PLASMA_FOUND)
    find_package_handle_standard_args(PLASMA DEFAULT_MSG
                                      PLASMA_LIBRARIES)
else()
    find_package_handle_standard_args(PLASMA DEFAULT_MSG
                                      PLASMA_LIBRARIES
                                      PLASMA_WORKS)
endif()
