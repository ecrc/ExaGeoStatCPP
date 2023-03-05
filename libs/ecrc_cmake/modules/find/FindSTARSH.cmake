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
# - Find STARSH include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(STARSH
#               [REQUIRED]             # Fail with error if starsh is not found
#               [COMPONENTS <comp1> <comp2> ...] # dependencies
#              )
#
#  STARSH depends on the following libraries:
#   - LAPACK
#   - LAPACKE
#   - BLAS
#   - CBLAS
#
#  COMPONENTS are optional libraries STARSH could be linked with,
#  Use it to drive detection of a specific compilation chain
#  COMPONENTS can be some of the following:
#   - no components are available for now: maybe STARSH in the future?
#
# Results are reported in variables:
#  STARSH_FOUND            - True if headers and requested libraries were found
#  STARSH_LINKER_FLAGS     - list of required linker flags (excluding -l and -L)
#  STARSH_INCLUDE_DIRS     - starsh include directories
#  STARSH_LIBRARY_DIRS     - Link directories for starsh libraries
#  STARSH_LIBRARIES        - starsh libraries
#  STARSH_INCLUDE_DIRS_DEP - starsh + dependencies include directories
#  STARSH_LIBRARY_DIRS_DEP - starsh + dependencies link directories
#  STARSH_LIBRARIES_DEP    - starsh libraries + dependencies
#
# The user can give specific paths where to find the libraries adding cmake
# options at configure (ex: cmake path/to/project -DSTARSH_DIR=path/to/starsh):
#  STARSH_DIR              - Where to find the base directory of starsh
#  STARSH_INCDIR           - Where to find the header files
#  STARSH_LIBDIR           - Where to find the library files
# The module can also look for the following environment variables if paths
# are not given as cmake variable: STARSH_DIR, STARSH_INCDIR, STARSH_LIBDIR
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


if(NOT STARSH_FOUND)
    set(STARSH_DIR "" CACHE PATH "Installation directory of STARSH library")
    if (NOT STARSH_FIND_QUIETLY)
        message(STATUS "A cache variable, namely STARSH_DIR, has been set to specify the install directory of STARSH")
    endif()
endif(NOT STARSH_FOUND)

# STARSH depends on LAPACKE anyway, try to find it
if (NOT LAPACKE_FOUND)
    if(STARSH_FIND_REQUIRED)
        find_package(LAPACKE REQUIRED)
    else()
        find_package(LAPACKE)
    endif()
endif()
# STARSH depends on CBLAS anyway, try to find it
if (NOT CBLAS_FOUND)
    if(STARSH_FIND_REQUIRED)
        find_package(CBLAS REQUIRED)
    else()
        find_package(CBLAS)
    endif()
endif()
# BLAS and LAPACK are searched by CBLAS and LAPACKE


set(ENV_STARSH_DIR "$ENV{STARSH_DIR}")
set(ENV_STARSH_INCDIR "$ENV{STARSH_INCDIR}")
set(ENV_STARSH_LIBDIR "$ENV{STARSH_LIBDIR}")
set(STARSH_GIVEN_BY_USER "FALSE")
if ( STARSH_DIR OR ( STARSH_INCDIR AND STARSH_LIBDIR) OR ENV_STARSH_DIR OR (ENV_STARSH_INCDIR AND ENV_STARSH_LIBDIR) )
    set(STARSH_GIVEN_BY_USER "TRUE")
endif()

# Optionally use pkg-config to detect include/library dirs (if pkg-config is available)
# -------------------------------------------------------------------------------------
include(FindPkgConfig)
find_package(PkgConfig QUIET)
if(PKG_CONFIG_EXECUTABLE AND NOT STARSH_GIVEN_BY_USER)

    pkg_search_module(STARSH starsh)
    if (NOT STARSH_FIND_QUIETLY)
        if (STARSH_FOUND AND STARSH_LIBRARIES)
            message(STATUS "Looking for STARSH - found using PkgConfig")
            #if(NOT STARSH_INCLUDE_DIRS)
            #    message("${Magenta}STARSH_INCLUDE_DIRS is empty using PkgConfig."
            #        "Perhaps the path to starsh headers is already present in your"
            #        "C(PLUS)_INCLUDE_PATH environment variable.${ColourReset}")
            #endif()
        else()
            message("${Magenta}Looking for STARSH - not found using PkgConfig. "
                "Perhaps you should add the directory containing starsh.pc "
                "to the PKG_CONFIG_PATH environment variable.${ColourReset}")
        endif()
    endif()

    if (STARSH_FIND_VERSION_EXACT)
        if( NOT (STARSH_FIND_VERSION_MAJOR STREQUAL STARSH_VERSION_MAJOR) OR
            NOT (STARSH_FIND_VERSION_MINOR STREQUAL STARSH_VERSION_MINOR) )
            if(NOT STARSH_FIND_QUIETLY)
                message(FATAL_ERROR
                        "STARSH version found is ${STARSH_VERSION_STRING} "
                        "when required is ${STARSH_FIND_VERSION}")
            endif()
        endif()
    else()
        # if the version found is older than the required then error
        if( (STARSH_FIND_VERSION_MAJOR STRGREATER STARSH_VERSION_MAJOR) OR
            (STARSH_FIND_VERSION_MINOR STRGREATER STARSH_VERSION_MINOR) )
            if(NOT STARSH_FIND_QUIETLY)
                message(FATAL_ERROR
                        "STARSH version found is ${STARSH_VERSION_STRING} "
                        "when required is ${STARSH_FIND_VERSION} or newer")
            endif()
        endif()
    endif()

    # if pkg-config is used: these variables are empty
    # the pkg_search_module call will set the following:
    # STARSH_LDFLAGS: all required linker flags
    # STARSH_CFLAGS:  all required cflags
    set(STARSH_INCLUDE_DIRS_DEP "")
    set(STARSH_LIBRARY_DIRS_DEP "")
    set(STARSH_LIBRARIES_DEP "")
    # replace it anyway: we should update it with dependencies given by pkg-config
    set(STARSH_INCLUDE_DIRS_DEP "${STARSH_INCLUDE_DIRS}")
    set(STARSH_LIBRARY_DIRS_DEP "${STARSH_LIBRARY_DIRS}")
    set(STARSH_LIBRARIES_DEP "${STARSH_LIBRARIES}")

endif(PKG_CONFIG_EXECUTABLE AND NOT STARSH_GIVEN_BY_USER)

# if STARSH is not found using pkg-config
if( (NOT PKG_CONFIG_EXECUTABLE) OR (PKG_CONFIG_EXECUTABLE AND NOT STARSH_FOUND) OR (STARSH_GIVEN_BY_USER) )

    if (NOT STARSH_FIND_QUIETLY)
        message(STATUS "Looking for STARSH - PkgConfig not used")
    endif()

    # Looking for include
    # -------------------

    # Add system include paths to search include
    # ------------------------------------------
    unset(_inc_env)
    set(ENV_STARSH_DIR "$ENV{STARSH_DIR}")
    set(ENV_STARSH_INCDIR "$ENV{STARSH_INCDIR}")
    if(ENV_STARSH_INCDIR)
        list(APPEND _inc_env "${ENV_STARSH_INCDIR}")
    elseif(ENV_STARSH_DIR)
        list(APPEND _inc_env "${ENV_STARSH_DIR}")
        list(APPEND _inc_env "${ENV_STARSH_DIR}/include")
        list(APPEND _inc_env "${ENV_STARSH_DIR}/include/starsh")
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


    # Try to find the starsh header in the given paths
    # -------------------------------------------------
    # call cmake macro to find the header path
    if(STARSH_INCDIR)
        set(STARSH_starsh.h_DIRS "STARSH_starsh.h_DIRS-NOTFOUND")
        find_path(STARSH_starsh.h_DIRS
          NAMES starsh.h
          HINTS ${STARSH_INCDIR})
    else()
        if(STARSH_DIR)
            set(STARSH_starsh.h_DIRS "STARSH_starsh.h_DIRS-NOTFOUND")
            find_path(STARSH_starsh.h_DIRS
              NAMES starsh.h
              HINTS ${STARSH_DIR}
              PATH_SUFFIXES "include" "include/starsh")
        else()
            set(STARSH_starsh.h_DIRS "STARSH_starsh.h_DIRS-NOTFOUND")
            find_path(STARSH_starsh.h_DIRS
              NAMES starsh.h
              HINTS ${_inc_env})
        endif()
    endif()
    mark_as_advanced(STARSH_starsh.h_DIRS)

    # If found, add path to cmake variable
    # ------------------------------------
    if (STARSH_starsh.h_DIRS)
        set(STARSH_INCLUDE_DIRS "${STARSH_starsh.h_DIRS}")
    else ()
        set(STARSH_INCLUDE_DIRS "STARSH_INCLUDE_DIRS-NOTFOUND")
        if(NOT STARSH_FIND_QUIETLY)
            message(STATUS "Looking for starsh -- starsh.h not found")
        endif()
    endif()


    # Looking for lib
    # ---------------

    # Add system library paths to search lib
    # --------------------------------------
    unset(_lib_env)
    set(ENV_STARSH_LIBDIR "$ENV{STARSH_LIBDIR}")
    if(ENV_STARSH_LIBDIR)
        list(APPEND _lib_env "${ENV_STARSH_LIBDIR}")
    elseif(ENV_STARSH_DIR)
        list(APPEND _lib_env "${ENV_STARSH_DIR}")
        list(APPEND _lib_env "${ENV_STARSH_DIR}/lib")
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

    # Try to find the starsh lib in the given paths
    # ----------------------------------------------

    # call cmake macro to find the lib path
    if(STARSH_LIBDIR)
        set(STARSH_starsh_LIBRARY "STARSH_starsh_LIBRARY-NOTFOUND")
        find_library(STARSH_starsh_LIBRARY
            NAMES starsh
            HINTS ${STARSH_LIBDIR})
    else()
        if(STARSH_DIR)
            set(STARSH_starsh_LIBRARY "STARSH_starsh_LIBRARY-NOTFOUND")
            find_library(STARSH_starsh_LIBRARY
                NAMES starsh
                HINTS ${STARSH_DIR}
                PATH_SUFFIXES lib lib32 lib64)
        else()
            set(STARSH_starsh_LIBRARY "STARSH_starsh_LIBRARY-NOTFOUND")
            find_library(STARSH_starsh_LIBRARY
                NAMES starsh
                HINTS ${_lib_env})
        endif()
    endif()
    mark_as_advanced(STARSH_starsh_LIBRARY)

    # If found, add path to cmake variable
    # ------------------------------------
    if (STARSH_starsh_LIBRARY)
        get_filename_component(starsh_lib_path "${STARSH_starsh_LIBRARY}" PATH)
        # set cmake variables
        set(STARSH_LIBRARIES    "${STARSH_starsh_LIBRARY}")
        set(STARSH_LIBRARY_DIRS "${starsh_lib_path}")
    else ()
        set(STARSH_LIBRARIES    "STARSH_LIBRARIES-NOTFOUND")
        set(STARSH_LIBRARY_DIRS "STARSH_LIBRARY_DIRS-NOTFOUND")
        if(NOT STARSH_FIND_QUIETLY)
            message(STATUS "Looking for starsh -- lib starsh not found")
        endif()
    endif ()

    # check a function to validate the find
    if (STARSH_LIBRARIES)

        set(REQUIRED_LDFLAGS)
        set(REQUIRED_INCDIRS)
        set(REQUIRED_LIBDIRS)
        set(REQUIRED_LIBS)

        # STARSH
        if (STARSH_INCLUDE_DIRS)
            set(REQUIRED_INCDIRS "${STARSH_INCLUDE_DIRS}")
        endif()
        if (STARSH_LIBRARY_DIRS)
            set(REQUIRED_LIBDIRS "${STARSH_LIBRARY_DIRS}")
        endif()
        set(REQUIRED_LIBS "${STARSH_LIBRARIES}")
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
        unset(STARSH_WORKS CACHE)
        include(CheckFunctionExists)
        check_function_exists(starsh_dgetrf STARSH_WORKS)
        mark_as_advanced(STARSH_WORKS)

        if(STARSH_WORKS)
            # save link with dependencies
            set(STARSH_LIBRARIES_DEP    "${REQUIRED_LIBS}")
            set(STARSH_LIBRARY_DIRS_DEP "${REQUIRED_LIBDIRS}")
            set(STARSH_INCLUDE_DIRS_DEP "${REQUIRED_INCDIRS}")
            set(STARSH_LINKER_FLAGS     "${REQUIRED_LDFLAGS}")
            list(REMOVE_DUPLICATES STARSH_LIBRARY_DIRS_DEP)
            list(REMOVE_DUPLICATES STARSH_INCLUDE_DIRS_DEP)
            list(REMOVE_DUPLICATES STARSH_LINKER_FLAGS)
        else()
            if(NOT STARSH_FIND_QUIETLY)
                message(STATUS "Looking for starsh : test of starsh_dgetrf with
                starsh, cblas, and lapack libraries fails")
                message(STATUS "CMAKE_REQUIRED_LIBRARIES: ${CMAKE_REQUIRED_LIBRARIES}")
                message(STATUS "CMAKE_REQUIRED_INCLUDES: ${CMAKE_REQUIRED_INCLUDES}")
                message(STATUS "Check in CMakeFiles/CMakeError.log to figure out why it fails")
            endif()
        endif()
        set(CMAKE_REQUIRED_INCLUDES)
        set(CMAKE_REQUIRED_FLAGS)
        set(CMAKE_REQUIRED_LIBRARIES)
    endif(STARSH_LIBRARIES)

endif( (NOT PKG_CONFIG_EXECUTABLE) OR (PKG_CONFIG_EXECUTABLE AND NOT STARSH_FOUND) OR (STARSH_GIVEN_BY_USER) )

if (STARSH_LIBRARIES)
    if (STARSH_LIBRARY_DIRS)
        set( first_lib_path "" )
        foreach(dir ${STARSH_LIBRARY_DIRS})
            if ("${dir}" MATCHES "starsh")
                set(first_lib_path "${dir}")
            endif()
        endforeach()
        if( NOT first_lib_path )
            list(GET STARSH_LIBRARY_DIRS 0 first_lib_path)
        endif()
    else()
        list(GET STARSH_LIBRARIES 0 first_lib)
        get_filename_component(first_lib_path "${first_lib}" PATH)
    endif()
    if (${first_lib_path} MATCHES "/lib(32|64)?$")
        string(REGEX REPLACE "/lib(32|64)?$" "" not_cached_dir "${first_lib_path}")
        set(STARSH_DIR_FOUND "${not_cached_dir}" CACHE PATH "Installation directory of STARSH library" FORCE)
    else()
        set(STARSH_DIR_FOUND "${first_lib_path}" CACHE PATH "Installation directory of STARSH library" FORCE)
    endif()
endif()

# check that STARSH has been found
# -------------------------------
include(FindPackageHandleStandardArgs)
if (PKG_CONFIG_EXECUTABLE AND STARSH_FOUND)
    find_package_handle_standard_args(STARSH DEFAULT_MSG
                                      STARSH_LIBRARIES)
else()
    find_package_handle_standard_args(STARSH DEFAULT_MSG
                                      STARSH_LIBRARIES
                                      STARSH_WORKS)
endif()
