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
#
###
#
# - Find HCORE include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(HCORE
#               [REQUIRED]             # Fail with error if HCORE is not found
#               [COMPONENTS <comp1> <comp2> ...] # dependencies
#              )
#
#  HCORE depends on the following libraries:
#   - CBLAS
#   - LAPACKE
#
#  COMPONENTS are optional libraries HCORE could be linked with,
#  Use it to drive detection of a specific compilation chain
#  COMPONENTS can be some of the following:
#
# This module finds headers and HCORE library.
# Results are reported in variables:
#  HCORE_FOUND            - True if headers and requested libraries were found
#  HCORE_LINKER_FLAGS     - list of required linker flags (excluding -l and -L)
#  HCORE_INCLUDE_DIRS     - HCORE include directories
#  HCORE_LIBRARY_DIRS     - Link directories for HCORE libraries
#  HCORE_INCLUDE_DIRS_DEP - HCORE + dependencies include directories
#  HCORE_LIBRARY_DIRS_DEP - HCORE + dependencies link directories
#  HCORE_LIBRARIES_DEP    - HCORE libraries + dependencies
# The user can give specific paths where to find the libraries adding cmake
# options at configure (ex: cmake path/to/project -DHCORE_DIR=path/to/HCORE):
#  HCORE_DIR              - Where to find the base directory of HCORE
#  HCORE_INCDIR           - Where to find the header files
#  HCORE_LIBDIR           - Where to find the library files
# The module can also look for the following environment variables if paths
# are not given as cmake variable: HCORE_DIR, HCORE_INCDIR, HCORE_LIBDIR

#=============================================================================
# Copyright 2012-2013 Inria
# Copyright 2012-2013 Emmanuel Agullo
# Copyright 2012-2013 Mathieu Faverge
# Copyright 2012      Cedric Castagnede
# Copyright 2013-2016 Florent Pruvost
# Copyright 2016-2019 Eduardo Gonzalez Fisher
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file ECRC-Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of ECRC, substitute the full
#  License text for the above reference.)


if (NOT HCORE_FOUND)
    set(HCORE_DIR "" CACHE PATH "Installation directory of HCORE library")
    if (NOT HCORE_FIND_QUIETLY)
        message(STATUS "A cache variable, namely HCORE_DIR, has been set to specify the install directory of HCORE")
    endif()
endif()

set(ENV_HCORE_DIR "$ENV{HCORE_DIR}")
set(ENV_HCORE_INCDIR "$ENV{HCORE_INCDIR}")
set(ENV_HCORE_LIBDIR "$ENV{HCORE_LIBDIR}")
set(HCORE_GIVEN_BY_USER "FALSE")
if ( HCORE_DIR OR ( HCORE_INCDIR AND HCORE_LIBDIR) OR ENV_HCORE_DIR OR (ENV_HCORE_INCDIR AND ENV_HCORE_LIBDIR) )
    set(HCORE_GIVEN_BY_USER "TRUE")
endif()

# Optionally use pkg-config to detect include/library dirs (if pkg-config is available)
# -------------------------------------------------------------------------------------
include(FindPkgConfig)
find_package(PkgConfig QUIET)
if(PKG_CONFIG_EXECUTABLE AND NOT HCORE_GIVEN_BY_USER)

    pkg_search_module(HCORE hcore)
    if (NOT HCORE_FIND_QUIETLY)
        if (HCORE_FOUND AND HCORE_LIBRARIES)
            message(STATUS "Looking for HCORE - found using PkgConfig")
            #if(NOT HCORE_INCLUDE_DIRS)
            # message("${Magenta}HCORE_INCLUDE_DIRS is empty using PkgConfig."
            # "Perhaps the path to HCORE headers is already present in your"
            # "C(PLUS)_INCLUDE_PATH environment variable.${ColourReset}")
            #endif()
        else()
            message(STATUS "${Magenta}Looking for HCORE - not found using PkgConfig."
            "\n     Perhaps you should add the directory containing hcore.pc"
            "\n     to the PKG_CONFIG_PATH environment variable.${ColourReset}")
        endif()
    endif()

    if (HCORE_FIND_VERSION_EXACT)
        if( NOT (HCORE_FIND_VERSION_MAJOR STREQUAL HCORE_VERSION_MAJOR) OR
            NOT (HCORE_FIND_VERSION_MINOR STREQUAL HCORE_VERSION_MINOR) )
            if(NOT HCORE_FIND_QUIETLY)
                message(FATAL_ERROR
                "HCORE version found is ${HCORE_VERSION_STRING}"
                "when required is ${HCORE_FIND_VERSION}")
            endif()
        endif()
    else()
        # if the version found is older than the required then error
        if( (HCORE_FIND_VERSION_MAJOR STRGREATER HCORE_VERSION_MAJOR) OR
            (HCORE_FIND_VERSION_MINOR STRGREATER HCORE_VERSION_MINOR) )
            if(NOT HCORE_FIND_QUIETLY)
                message(FATAL_ERROR
                "HCORE version found is ${HCORE_VERSION_STRING}"
                "when required is ${HCORE_FIND_VERSION} or newer")
            endif()
        endif()
    endif()

    set(HCORE_INCLUDE_DIRS_DEP "${HCORE_INCLUDE_DIRS}")
    set(HCORE_LIBRARY_DIRS_DEP "${HCORE_LIBRARY_DIRS}")
    set(HCORE_LIBRARIES_DEP "${HCORE_LIBRARIES}")

endif(PKG_CONFIG_EXECUTABLE AND NOT HCORE_GIVEN_BY_USER)

if( (NOT PKG_CONFIG_EXECUTABLE) OR (PKG_CONFIG_EXECUTABLE AND NOT HCORE_FOUND) OR (HCORE_GIVEN_BY_USER) )

    if (NOT HCORE_FIND_QUIETLY)
        message(STATUS "Looking for HCORE - PkgConfig not used")
    endif()

    # Dependencies detection
    # ----------------------

    if (NOT HCORE_FIND_QUIETLY)
        message(STATUS "Looking for HCORE - Try to detect pthread")
    endif()
    if (HCORE_FIND_REQUIRED)
        find_package(Threads REQUIRED)
    else()
        find_package(Threads)
    endif()
    set(HCORE_EXTRA_LIBRARIES "")
    if( THREADS_FOUND )
        list(APPEND HCORE_EXTRA_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})
    endif ()

    # Add math library to the list of extra
    # it normally exists on all common systems provided with a C compiler
    if (NOT HCORE_FIND_QUIETLY)
        message(STATUS "Looking for HCORE - Try to detect libm")
    endif()
    set(HCORE_M_LIBRARIES "")
    if(UNIX OR WIN32)
        find_library(
            HCORE_M_m_LIBRARY
            NAMES m
            )
        mark_as_advanced(HCORE_M_m_LIBRARY)
        if (HCORE_M_m_LIBRARY)
            list(APPEND HCORE_M_LIBRARIES "${HCORE_M_m_LIBRARY}")
            list(APPEND HCORE_EXTRA_LIBRARIES "${HCORE_M_m_LIBRARY}")
        else()
            if (HCORE_FIND_REQUIRED)
                message(FATAL_ERROR "Could NOT find libm on your system."
                "Are you sure to a have a C compiler installed?")
            endif()
        endif()
    endif()

    # Try to find librt (libposix4 - POSIX.1b Realtime Extensions library)
    # on Unix systems except Apple ones because it does not exist on it
    if (NOT HCORE_FIND_QUIETLY)
        message(STATUS "Looking for HCORE - Try to detect librt")
    endif()
    set(HCORE_RT_LIBRARIES "")
    if(UNIX AND NOT APPLE)
        find_library(
            HCORE_RT_rt_LIBRARY
            NAMES rt
            )
        mark_as_advanced(HCORE_RT_rt_LIBRARY)
        if (HCORE_RT_rt_LIBRARY)
            list(APPEND HCORE_RT_LIBRARIES "${HCORE_RT_rt_LIBRARY}")
            list(APPEND HCORE_EXTRA_LIBRARIES "${HCORE_RT_rt_LIBRARY}")
        else()
            if (HCORE_FIND_REQUIRED)
                message(FATAL_ERROR "Could NOT find librt on your system")
            endif()
        endif()
    endif()

    # HCORE depends on CBLAS
    #---------------------------
    if (NOT HCORE_FIND_QUIETLY)
        message(STATUS "Looking for HCORE - Try to detect CBLAS (depends on BLAS)")
    endif()
    if (HCORE_FIND_REQUIRED)
        find_package(CBLAS REQUIRED)
    else()
        find_package(CBLAS)
    endif()

    # HCORE depends on LAPACKE
    #-----------------------------

    # standalone version of lapacke seems useless for now
    # let the comment in case we meet some problems of non existing lapacke
    # functions in lapack library such as mkl, acml, ...
    #set(LAPACKE_STANDALONE TRUE)
    if (NOT HCORE_FIND_QUIETLY)
        message(STATUS "Looking for HCORE - Try to detect LAPACKE (depends on LAPACK)")
    endif()
    if (HCORE_FIND_REQUIRED)
        find_package(LAPACKE REQUIRED)
    else()
        find_package(LAPACKE)
    endif()

    # Looking for include
    # -------------------

    # Add system include paths to search include
    # ------------------------------------------
    unset(_inc_env)
    set(ENV_HCORE_DIR "$ENV{HCORE_DIR}")
    set(ENV_HCORE_INCDIR "$ENV{HCORE_INCDIR}")
    if(ENV_HCORE_INCDIR)
        list(APPEND _inc_env "${ENV_HCORE_INCDIR}")
    elseif(ENV_HCORE_DIR)
        list(APPEND _inc_env "${ENV_HCORE_DIR}")
        list(APPEND _inc_env "${ENV_HCORE_DIR}/include")
        list(APPEND _inc_env "${ENV_HCORE_DIR}/include/hcore")
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


    # Try to find the hcore header in the given paths
    # ---------------------------------------------------
    # call cmake macro to find the header path
    if(HCORE_INCDIR)
        set(hcore.h_DIRS "hcore.h_DIRS-NOTFOUND")
        find_path(hcore.h_DIRS
            NAMES hcore.h
            HINTS ${HCORE_INCDIR})
    else()
        if(HCORE_DIR)
            set(hcore.h_DIRS "hcore.h_DIRS-NOTFOUND")
            find_path(hcore.h_DIRS
            NAMES hcore.h
            HINTS ${HCORE_DIR}
            PATH_SUFFIXES "include" "include/hcore")
        else()
            set(hcore.h_DIRS "hcore.h_DIRS-NOTFOUND")
            find_path(hcore.h_DIRS
            NAMES hcore.h
            HINTS ${_inc_env}
            PATH_SUFFIXES "hcore")
        endif()
    endif()
    mark_as_advanced(hcore.h_DIRS)

    # If found, add path to cmake variable
    # ------------------------------------
    if (hcore.h_DIRS)
        set(HCORE_INCLUDE_DIRS "${hcore.h_DIRS}")
    else ()
        set(HCORE_INCLUDE_DIRS "HCORE_INCLUDE_DIRS-NOTFOUND")
        if(NOT HCORE_FIND_QUIETLY)
            message(STATUS "Looking for HCORE -- hcore.h not found")
        endif()
    endif()


    # Looking for lib
    # ---------------

    # Add system library paths to search lib
    # --------------------------------------
    unset(_lib_env)
    set(ENV_HCORE_LIBDIR "$ENV{HCORE_LIBDIR}")
    if(ENV_HCORE_LIBDIR)
        list(APPEND _lib_env "${ENV_HCORE_LIBDIR}")
    elseif(ENV_HCORE_DIR)
        list(APPEND _lib_env "${ENV_HCORE_DIR}")
        list(APPEND _lib_env "${ENV_HCORE_DIR}/lib")
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

    # Try to find the HCORE lib in the given paths
    # ------------------------------------------------

    # create list of libs to find
    set(HCORE_libs_to_find "hcore")
    list(APPEND HCORE_libs_to_find "coreblas")

    # call cmake macro to find the lib path
    if(HCORE_LIBDIR)
        foreach(hcore_lib ${HCORE_libs_to_find})
            set(HCORE_${hcore_lib}_LIBRARY "HCORE_${hcore_lib}_LIBRARY-NOTFOUND")
            find_library(HCORE_${hcore_lib}_LIBRARY
            NAMES ${hcore_lib}
            HINTS ${HCORE_LIBDIR})
        endforeach()
    else()
        if(HCORE_DIR)
            foreach(hcore_lib ${HCORE_libs_to_find})
                set(HCORE_${hcore_lib}_LIBRARY "HCORE_${hcore_lib}_LIBRARY-NOTFOUND")
                find_library(HCORE_${hcore_lib}_LIBRARY
                NAMES ${hcore_lib}
                HINTS ${HCORE_DIR}
                PATH_SUFFIXES lib lib32 lib64)
            endforeach()
        else()
            foreach(hcore_lib ${HCORE_libs_to_find})
                set(HCORE_${hcore_lib}_LIBRARY "HCORE_${hcore_lib}_LIBRARY-NOTFOUND")
                find_library(HCORE_${hcore_lib}_LIBRARY
                NAMES ${hcore_lib}
                HINTS ${_lib_env})
            endforeach()
        endif()
    endif()

    # If found, add path to cmake variable
    # ------------------------------------
    foreach(hcore_lib ${HCORE_libs_to_find})

        get_filename_component(${hcore_lib}_lib_path ${HCORE_${hcore_lib}_LIBRARY} PATH)
        # set cmake variables (respects naming convention)
        if (HCORE_LIBRARIES)
            list(APPEND HCORE_LIBRARIES "${HCORE_${hcore_lib}_LIBRARY}")
        else()
            set(HCORE_LIBRARIES "${HCORE_${hcore_lib}_LIBRARY}")
        endif()
        if (HCORE_LIBRARY_DIRS)
            list(APPEND HCORE_LIBRARY_DIRS "${${hcore_lib}_lib_path}")
        else()
            set(HCORE_LIBRARY_DIRS "${${hcore_lib}_lib_path}")
        endif()
        mark_as_advanced(HCORE_${hcore_lib}_LIBRARY)

    endforeach(hcore_lib ${HCORE_libs_to_find})

    # check a function to validate the find
    if(HCORE_LIBRARIES)

        set(REQUIRED_LDFLAGS)
        set(REQUIRED_INCDIRS)
        set(REQUIRED_LIBDIRS)
        set(REQUIRED_LIBS)

        # HCORE
        if (HCORE_INCLUDE_DIRS)
            set(REQUIRED_INCDIRS "${HCORE_INCLUDE_DIRS}")
        endif()
        foreach(libdir ${HCORE_LIBRARY_DIRS})
            if (libdir)
                list(APPEND REQUIRED_LIBDIRS "${libdir}")
            endif()
        endforeach()
        set(REQUIRED_LIBS "${HCORE_LIBRARIES}")
        # LAPACKE
        if (LAPACKE_FOUND)
            if (LAPACKE_INCLUDE_DIRS_DEP)
                list(APPEND REQUIRED_INCDIRS "${LAPACKE_INCLUDE_DIRS_DEP}")
            elseif (LAPACKE_INCLUDE_DIRS)
                list(APPEND REQUIRED_INCDIRS "${LAPACKE_INCLUDE_DIRS}")
            endif()
            if(LAPACKE_LIBRARY_DIRS_DEP)
                list(APPEND REQUIRED_LIBDIRS "${LAPACKE_LIBRARY_DIRS_DEP}")
            elseif(LAPACKE_LIBRARY_DIRS)
                list(APPEND REQUIRED_LIBDIRS "${LAPACKE_LIBRARY_DIRS}")
            endif()
            if (LAPACKE_LIBRARIES_DEP)
                list(APPEND REQUIRED_LIBS "${LAPACKE_LIBRARIES_DEP}")
            elseif(LAPACKE_LIBRARIES)
                list(APPEND REQUIRED_LIBS "${LAPACKE_LIBRARIES}")
            endif()
            if (LAPACK_LINKER_FLAGS)
                list(APPEND REQUIRED_LDFLAGS "${LAPACK_LINKER_FLAGS}")
            endif()
        endif()
        # CBLAS
        if (CBLAS_FOUND)
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
        endif()
        # EXTRA LIBS such that pthread, m, rt
        list(APPEND REQUIRED_LIBS ${HCORE_EXTRA_LIBRARIES})

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
        unset(HCORE_WORKS CACHE)
        include(CheckFunctionExists)
        check_function_exists(HCORE_Init HCORE_WORKS)
        mark_as_advanced(HCORE_WORKS)

        if(HCORE_WORKS)
            # save link with dependencies
            set(HCORE_LIBRARIES_DEP "${REQUIRED_LIBS}")
            set(HCORE_LIBRARY_DIRS_DEP "${REQUIRED_LIBDIRS}")
            set(HCORE_INCLUDE_DIRS_DEP "${REQUIRED_INCDIRS}")
            set(HCORE_LINKER_FLAGS "${REQUIRED_LDFLAGS}")
            list(REMOVE_DUPLICATES HCORE_LIBRARY_DIRS_DEP)
            list(REMOVE_DUPLICATES HCORE_INCLUDE_DIRS_DEP)
            list(REMOVE_DUPLICATES HCORE_LINKER_FLAGS)
        else()
            if(NOT HCORE_FIND_QUIETLY)
                message(STATUS "Looking for HCORE : test of HCORE_Init fails")
                message(STATUS "CMAKE_REQUIRED_LIBRARIES: ${CMAKE_REQUIRED_LIBRARIES}")
                message(STATUS "CMAKE_REQUIRED_INCLUDES: ${CMAKE_REQUIRED_INCLUDES}")
                message(STATUS "Check in CMakeFiles/CMakeError.log to figure out why it fails")
                message(STATUS "Maybe HCORE is linked with specific libraries. "
                "See the explanation in FindHCORE.cmake.")
            endif()
        endif()
        set(CMAKE_REQUIRED_INCLUDES)
        set(CMAKE_REQUIRED_FLAGS)
        set(CMAKE_REQUIRED_LIBRARIES)
    endif(HCORE_LIBRARIES)

endif( (NOT PKG_CONFIG_EXECUTABLE) OR (PKG_CONFIG_EXECUTABLE AND NOT HCORE_FOUND) OR (HCORE_GIVEN_BY_USER) )

if (HCORE_LIBRARIES)
    if (HCORE_LIBRARY_DIRS)
        set( first_lib_path "" )
        foreach(dir ${HCORE_LIBRARY_DIRS})
            if ("${dir}" MATCHES "hcore")
                set(first_lib_path "${dir}")
            endif()
        endforeach()
        if( NOT first_lib_path )
            list(GET HCORE_LIBRARY_DIRS 0 first_lib_path)
        endif()
    else()
        list(GET HCORE_LIBRARIES 0 first_lib)
        get_filename_component(first_lib_path "${first_lib}" PATH)
    endif()
    if (${first_lib_path} MATCHES "/lib(32|64)?$")
        string(REGEX REPLACE "/lib(32|64)?$" "" not_cached_dir "${first_lib_path}")
        set(HCORE_DIR_FOUND "${not_cached_dir}" CACHE PATH "Installation directory of HCORE library" FORCE)
    else()
        set(HCORE_DIR_FOUND "${first_lib_path}" CACHE PATH "Installation directory of HCORE library" FORCE)
    endif()
endif()
mark_as_advanced(HCORE_DIR)
mark_as_advanced(HCORE_DIR_FOUND)

# check that HCORE has been found
# ---------------------------------
include(FindPackageHandleStandardArgs)
if (PKG_CONFIG_EXECUTABLE AND HCORE_FOUND)
    find_package_handle_standard_args(HCORE DEFAULT_MSG
        HCORE_LIBRARIES)
else()
    find_package_handle_standard_args(HCORE DEFAULT_MSG
        HCORE_LIBRARIES
        HCORE_WORKS)
endif()
