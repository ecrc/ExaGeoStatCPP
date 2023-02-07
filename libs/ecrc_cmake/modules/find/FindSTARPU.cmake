###
#
# @copyright (c) 2009-2014 The University of Tennessee and The University
#                          of Tennessee Research Foundation.
#                          All rights reserved.
# @copyright (c) 2012-2015 Inria. All rights reserved.
# @copyright (c) 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
# @copyright (c) 2022 King Abdullah University of Science and Technology (KAUST).
#                     All rights reserved.
#
###
#
# - Find STARPU include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(STARPU
#               [version] [EXACT]      # Minimum or EXACT version e.g. 1.1
#               [REQUIRED]             # Fail with error if starpu is not found
#               [COMPONENTS <comp1> <comp2> ...] # dependencies
#              )
#
#  COMPONENTS are optional libraries STARPU could be linked with,
#  Use it to drive detection of a specific compilation chain
#  COMPONENTS can be some of the following:
#   - HWLOC: to activate the detection of StarPU linked with HWLOC
#   - CUDA: to activate the detection of StarPU linked with CUDA
#   - MPI: to activate the detection of StarPU linked with MPI
#   - BLAS: to activate the detection of StarPU linked with BLAS
#   - MAGMA: to activate the detection of StarPU linked with MAGMA
#   - FXT: to activate the detection of StarPU linked with FxT
#   - SIMGRID: to activate the detection of StarPU linked with SimGrid
#
# Results are reported in variables:
#  STARPU_FOUND                  - True if headers and requested libraries were found
#  STARPU_C_FLAGS                - list of required compilation flags (excluding -I)
#  STARPU_LINKER_FLAGS           - list of required linker flags (excluding -l and -L)
#  STARPU_INCLUDE_DIRS           - starpu include directories
#  STARPU_LIBRARY_DIRS           - Link directories for starpu libraries
#  STARPU_LIBRARIES              - starpu libraries
#  STARPU_SHM_LIBRARIES          - starpu libraries without libstarpumpi
#  STARPU_MPI_LIBRARIES          - starpu libraries with libstarpumpi
#  STARPU_INCLUDE_DIRS_DEP       - starpu + dependencies include directories
#  STARPU_LIBRARY_DIRS_DEP       - starpu + dependencies link directories
#  STARPU_LIBRARIES_DEP          - starpu libraries + dependencies
#  STARPU_VERSION_STRING         - A human-readable string containing the version of the package found
#  STARPU_VERSION_MAJOR          - The major version of the package found
#  STARPU_VERSION_MINOR          - The minor version of the package found
# The user can give specific paths where to find the libraries adding cmake
# options at configure (ex: cmake path/to/project -DSTARPU=path/to/starpu):
#  STARPU_DIR                    - Where to find the base directory of starpu
#  STARPU_INCDIR                 - Where to find the header files
#  STARPU_LIBDIR                 - Where to find the library files
# The module can also look for the following environment variables if paths
# are not given as cmake variable: STARPU_DIR, STARPU_INCDIR, STARPU_LIBDIR

#=============================================================================
# Copyright 2012-2013 Inria
# Copyright 2012-2013 Emmanuel Agullo
# Copyright 2012-2013 Mathieu Faverge
# Copyright 2012      Cedric Castagnede
# Copyright 2013      Florent Pruvost
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file ECRC-Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of Ecrc, substitute the full
#  License text for the above reference.)

if (NOT STARPU_FOUND)
  set(STARPU_DIR "" CACHE PATH "Installation directory of STARPU library")
  if (NOT STARPU_FIND_QUIETLY)
    message(STATUS "A cache variable, namely STARPU_DIR, has been set to specify the install directory of STARPU")
  endif()
endif()

# STARPU may depend on other packages (HWLOC, MPI, CUDA, BLAS, MAGMA, FXT, ...)
# try to find them if specified as COMPONENTS during the call
set(STARPU_LOOK_FOR_HWLOC FALSE)
set(STARPU_LOOK_FOR_CUDA FALSE)
set(STARPU_LOOK_FOR_MPI FALSE)
set(STARPU_LOOK_FOR_BLAS FALSE)
set(STARPU_LOOK_FOR_MAGMA FALSE)
set(STARPU_LOOK_FOR_FXT FALSE)
set(STARPU_LOOK_FOR_SIMGRID FALSE)

if( STARPU_FIND_COMPONENTS )
  foreach( component ${STARPU_FIND_COMPONENTS} )
    if(${component} STREQUAL "HWLOC")
      set(STARPU_LOOK_FOR_HWLOC TRUE)
    elseif(${component} STREQUAL "CUDA")
      set(STARPU_LOOK_FOR_CUDA TRUE)
    elseif(${component} STREQUAL "MPI")
      set(STARPU_LOOK_FOR_MPI TRUE)
    elseif(${component} STREQUAL "BLAS")
      set(STARPU_LOOK_FOR_BLAS TRUE)
    elseif(${component} STREQUAL "MAGMA")
      set(STARPU_LOOK_FOR_MAGMA TRUE)
    elseif(${component} STREQUAL "FXT")
      set(STARPU_LOOK_FOR_FXT TRUE)
    elseif(${component} STREQUAL "SIMGRID")
      set(STARPU_LOOK_FOR_SIMGRID TRUE)
    endif()
  endforeach()
endif()

# STARPU may depend on pthread, try to find it
find_package(Threads)
if( THREADS_FOUND )
  list(APPEND STARPU_EXTRA_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})
endif ()
# STARPU may depend on libm, try to find it
find_library(M_m_LIBRARY NAMES m)
mark_as_advanced(M_m_LIBRARY)
if( M_m_LIBRARY )
  list(APPEND STARPU_EXTRA_LIBRARIES ${M_m_LIBRARY})
endif ()
# STARPU may depend on librt, try to find it
find_library(RT_rt_LIBRARY NAMES m)
mark_as_advanced(RT_rt_LIBRARY)
if( RT_rt_LIBRARY )
  list(APPEND STARPU_EXTRA_LIBRARIES ${RT_rt_LIBRARY})
endif ()

# STARPU may depend on HWLOC, try to find it
if (NOT HWLOC_FOUND AND STARPU_LOOK_FOR_HWLOC)
  if (STARPU_FIND_REQUIRED AND STARPU_FIND_REQUIRED_HWLOC)
    find_package(HWLOC REQUIRED)
  else()
    find_package(HWLOC)
  endif()
endif()

# STARPU may depend on CUDA, try to find it
if (NOT CUDA_FOUND AND STARPU_LOOK_FOR_CUDA)
  if (STARPU_FIND_REQUIRED AND STARPU_FIND_REQUIRED_CUDA)
    find_package(CUDA REQUIRED)
  else()
    find_package(CUDA)
  endif()
  if (CUDA_FOUND)
    mark_as_advanced(CUDA_BUILD_CUBIN)
    mark_as_advanced(CUDA_BUILD_EMULATION)
    mark_as_advanced(CUDA_SDK_ROOT_DIR)
    mark_as_advanced(CUDA_TOOLKIT_ROOT_DIR)
    mark_as_advanced(CUDA_VERBOSE_BUILD)
  endif()
endif()

# STARPU may depend on MPI, try to find it
if (NOT MPI_FOUND AND STARPU_LOOK_FOR_MPI)
  if (STARPU_FIND_REQUIRED AND STARPU_FIND_REQUIRED_MPI)
    find_package(MPI REQUIRED)
  else()
    find_package(MPI)
  endif()
  if (MPI_FOUND)
    mark_as_advanced(MPI_LIBRARY)
    mark_as_advanced(MPI_EXTRA_LIBRARY)
  endif()
endif()

# STARPU may depend on BLAS, try to find it
if (NOT BLAS_FOUND AND STARPU_LOOK_FOR_BLAS)
  if (STARPU_FIND_REQUIRED AND STARPU_FIND_REQUIRED_BLAS)
    find_package(BLASEXT REQUIRED)
  else()
    find_package(BLASEXT)
  endif()
endif()

# STARPU may depend on MAGMA, try to find it
if (NOT MAGMA_FOUND AND STARPU_LOOK_FOR_MAGMA)
  if (STARPU_FIND_REQUIRED AND STARPU_FIND_REQUIRED_MAGMA)
    find_package(MAGMA REQUIRED)
  else()
    find_package(MAGMA)
  endif()
endif()

# STARPU may depend on FXT, try to find it
if (NOT FXT_FOUND AND STARPU_LOOK_FOR_FXT)
  if (STARPU_FIND_REQUIRED AND STARPU_FIND_REQUIRED_FXT)
    find_package(FXT REQUIRED)
  else()
    find_package(FXT)
  endif()
endif()

# STARPU may depend on SIMGRID, try to find it
if (NOT SIMGRID_FOUND AND STARPU_LOOK_FOR_SIMGRID)
  if (STARPU_FIND_REQUIRED AND STARPU_FIND_REQUIRED_SIMGRID)
    find_package(SIMGRID REQUIRED)
  else()
    find_package(SIMGRID)
  endif()
endif()

set(ENV_STARPU_DIR "$ENV{STARPU_DIR}")
set(ENV_STARPU_INCDIR "$ENV{STARPU_INCDIR}")
set(ENV_STARPU_LIBDIR "$ENV{STARPU_LIBDIR}")
set(STARPU_GIVEN_BY_USER "FALSE")
if ( STARPU_DIR OR ( STARPU_INCDIR AND STARPU_LIBDIR) OR ENV_STARPU_DIR OR (ENV_STARPU_INCDIR AND ENV_STARPU_LIBDIR) )
  set(STARPU_GIVEN_BY_USER "TRUE")
endif()

# Optionally use pkg-config to detect include/library dirs (if pkg-config is available)
# -------------------------------------------------------------------------------------
include(FindPkgConfig)
find_package(PkgConfig QUIET)

if(PKG_CONFIG_EXECUTABLE AND NOT STARPU_GIVEN_BY_USER)

  pkg_search_module(STARPU_SHM libstarpu)
  set(STARPU_INCLUDE_DIRS "${STARPU_SHM_INCLUDE_DIRS}")
  set(STARPU_LIBRARY_DIRS "${STARPU_SHM_LIBRARY_DIRS}")
  if(STARPU_LOOK_FOR_MPI AND MPI_FOUND)
    pkg_search_module(STARPU_MPI libstarpumpi)
  endif()
  if (NOT STARPU_FIND_QUIETLY)
    if (STARPU_SHM_FOUND AND STARPU_SHM_LIBRARIES)
      message(STATUS "Looking for STARPU - found using PkgConfig")
      #if(NOT STARPU_SHM_INCLUDE_DIRS)
      #    message("${Magenta}STARPU_SHM_INCLUDE_DIRS is empty using PkgConfig."
      #        "Perhaps the path to starpu headers is already present in your"
      #        "C(PLUS)_INCLUDE_PATH environment variable.${ColourReset}")
      #endif()
      set(STARPU_VERSION_STRING "${STARPU_SHM_VERSION}")
      string(REPLACE "." ";" STARPU_VERSION_STRING_LIST ${STARPU_VERSION_STRING})
      list(GET STARPU_VERSION_STRING_LIST 0 STARPU_VERSION_MAJOR)
      list(GET STARPU_VERSION_STRING_LIST 1 STARPU_VERSION_MINOR)
    else()
      message(STATUS "${Magenta}Looking for STARPU - not found using PkgConfig."
	"\n   Perhaps you should add the directory containing libstarpu.pc"
	"\n   to the PKG_CONFIG_PATH environment variable.${ColourReset}")
    endif()
    if (STARPU_LOOK_FOR_MPI)
      if (STARPU_MPI_FOUND AND STARPU_MPI_LIBRARIES)
	message(STATUS "Looking for STARPU MPI - found using PkgConfig")
      else()
	message(STATUS "${Magenta}Looking for STARPU MPI - not found using PkgConfig."
	  " Perhaps you should add the directory containing libstarpumpi.pc to"
	  " the PKG_CONFIG_PATH environment variable.${ColourReset}")
      endif()
    endif()
  endif()

  if (STARPU_FIND_VERSION_EXACT)
    if( NOT (STARPU_FIND_VERSION_MAJOR STREQUAL STARPU_VERSION_MAJOR) OR
	NOT (STARPU_FIND_VERSION_MINOR STREQUAL STARPU_VERSION_MINOR) )
      if(STARPU_FIND_REQUIRED AND NOT STARPU_FIND_QUIETLY)
	message(FATAL_ERROR
	  "STARPU version found is ${STARPU_VERSION_STRING}"
	  " when required is ${STARPU_FIND_VERSION}")
      endif()
    endif()
  else()
    # if the version found is older than the required then error
    if( (STARPU_FIND_VERSION_MAJOR STRGREATER STARPU_VERSION_MAJOR) OR
	(STARPU_FIND_VERSION_MINOR STRGREATER STARPU_VERSION_MINOR) )
      if(STARPU_FIND_REQUIRED AND NOT STARPU_FIND_QUIETLY)
	message(FATAL_ERROR
	  "STARPU version found is ${STARPU_VERSION_STRING}"
	  " when required is ${STARPU_FIND_VERSION} or newer")
      endif()
    endif()
  endif()

  if(STARPU_MPI_LIBRARIES)
    if (STARPU_LOOK_FOR_SIMGRID)
      # Cmake does not fetch explicit libfxt.a static paths from pkg-config...
      find_package(FXT)
      string(REGEX MATCH "[^;]*/libfxt.a" FXT_STATIC_LIB "${STARPU_MPI_LDFLAGS_OTHER}")
      list(APPEND STARPU_MPI_LIBRARIES "${FXT_STATIC_LIB}")
    endif()
    set(STARPU_LIBRARIES "${STARPU_MPI_LIBRARIES}")
    set(STARPU_LINKER_FLAGS "${STARPU_MPI_LDFLAGS_OTHER}")
  elseif(STARPU_SHM_LIBRARIES)
    set(STARPU_LIBRARIES "${STARPU_SHM_LIBRARIES}")
    set(STARPU_LINKER_FLAGS "${STARPU_SHM_LDFLAGS_OTHER}")
  else()
    set(STARPU_LIBRARIES "STARPU_LIBRARIES-NOTFOUND")
  endif()
  set(STARPU_INCLUDE_DIRS_DEP "${STARPU_INCLUDE_DIRS}")
  set(STARPU_LIBRARY_DIRS_DEP "${STARPU_LIBRARY_DIRS}")
  set(STARPU_LIBRARIES_DEP "${STARPU_LIBRARIES}")
  if (STARPU_SHM_FOUND)
    set(STARPU_FOUND "TRUE")
  else()
    set(STARPU_FOUND "FALSE")
  endif()
  if (STARPU_LOOK_FOR_MPI AND NOT STARPU_MPI_FOUND)
    set(STARPU_FOUND "FALSE")
  endif()
endif(PKG_CONFIG_EXECUTABLE AND NOT STARPU_GIVEN_BY_USER)


if( (NOT PKG_CONFIG_EXECUTABLE) OR (PKG_CONFIG_EXECUTABLE AND NOT STARPU_FOUND) OR (STARPU_GIVEN_BY_USER) )

  # Looking for include
  # -------------------

  # Add system include paths to search include
  # ------------------------------------------
  unset(_inc_env)
  set(ENV_STARPU_DIR "$ENV{STARPU_DIR}")
  set(ENV_STARPU_INCDIR "$ENV{STARPU_INCDIR}")
  if(ENV_STARPU_INCDIR)
    list(APPEND _inc_env "${ENV_STARPU_INCDIR}")
  elseif(ENV_STARPU_DIR)
    list(APPEND _inc_env "${ENV_STARPU_DIR}")
    list(APPEND _inc_env "${ENV_STARPU_DIR}/include")
    list(APPEND _inc_env "${ENV_STARPU_DIR}/include/starpu/1.0")
    list(APPEND _inc_env "${ENV_STARPU_DIR}/include/starpu/1.1")
    list(APPEND _inc_env "${ENV_STARPU_DIR}/include/starpu/1.2")
    list(APPEND _inc_env "${ENV_STARPU_DIR}/include/starpu/1.3")
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

  # Try to find the version of StarPU in starpu_config.h file
  set(STARPU_hdrs_to_find "starpu_config.h")

  # call cmake macro to find the header path
  if(STARPU_INCDIR)
    foreach(starpu_hdr ${STARPU_hdrs_to_find})
      set(STARPU_${starpu_hdr}_INCLUDE_DIRS "STARPU_${starpu_hdr}_INCLUDE_DIRS-NOTFOUND")
      find_path(STARPU_${starpu_hdr}_INCLUDE_DIRS
	NAMES ${starpu_hdr}
	HINTS ${STARPU_INCDIR})
    endforeach()
  else()
    if(STARPU_DIR)
      foreach(starpu_hdr ${STARPU_hdrs_to_find})
	set(STARPU_${starpu_hdr}_INCLUDE_DIRS "STARPU_${starpu_hdr}_INCLUDE_DIRS-NOTFOUND")
	find_path(STARPU_${starpu_hdr}_INCLUDE_DIRS
	  NAMES ${starpu_hdr}
	  HINTS ${STARPU_DIR}
	  PATH_SUFFIXES "include"
	  "include/starpu/1.0"
	  "include/starpu/1.1"
	  "include/starpu/1.2"
	  "include/starpu/1.3")
      endforeach()
    else()
      foreach(starpu_hdr ${STARPU_hdrs_to_find})
	set(STARPU_${starpu_hdr}_INCLUDE_DIRS "STARPU_${starpu_hdr}_INCLUDE_DIRS-NOTFOUND")
	find_path(STARPU_${starpu_hdr}_INCLUDE_DIRS
	  NAMES ${starpu_hdr}
	  HINTS ${_inc_env}
	  PATH_SUFFIXES
	  "starpu/1.0"
	  "starpu/1.1"
	  "starpu/1.2"
	  "starpu/1.3")
      endforeach()
    endif()
  endif()
  mark_as_advanced(STARPU_starpu_config.h_INCLUDE_DIRS)

  ###
  #
  # GET_VERSION: Get the version of the software by parsing a file
  #
  ###
  MACRO(GET_VERSION _PACKAGE _filepath)

    #message(STATUS "Looking for ${_PACKAGE} version in the file ${_filepath}")
    file(READ "${_filepath}" _file)
    string(REGEX REPLACE
      "(.*)define([ \t]*)${_PACKAGE}_MAJOR_VERSION([ \t]*)([0-9]+)(.*)"
      "\\4" ${_PACKAGE}_VERSION_MAJOR "${_file}")
    string(REGEX REPLACE
      "(.*)define([ \t]*)${_PACKAGE}_MINOR_VERSION([ \t]*)([0-9]+)(.*)"
      "\\4" ${_PACKAGE}_VERSION_MINOR "${_file}")
    set(${_PACKAGE}_VERSION_STRING
      "${${_PACKAGE}_VERSION_MAJOR}.${${_PACKAGE}_VERSION_MINOR}" CACHE PATH "StarPU version found in starpu_config.h")
    #message(STATUS "${_PACKAGE}_VERSION_MAJOR = -${${_PACKAGE}_VERSION_MAJOR}-")
    #message(STATUS "${_PACKAGE}_VERSION_MINOR = -${${_PACKAGE}_VERSION_MINOR}-")

  ENDMACRO(GET_VERSION)

  # Find the version of StarPU in starpu_config.h file
  # remark: the version is defined in this file since the STARPU 1.0 version
  if (STARPU_starpu_config.h_INCLUDE_DIRS)
    GET_VERSION("STARPU" "${STARPU_starpu_config.h_INCLUDE_DIRS}/starpu_config.h")
    if (STARPU_VERSION_MAJOR AND STARPU_VERSION_MINOR)
      if (STARPU_FIND_VERSION_EXACT)
	if (STARPU_FIND_VERSION_MAJOR AND STARPU_FIND_VERSION_MINOR)
	  if( NOT (STARPU_FIND_VERSION_MAJOR STREQUAL STARPU_VERSION_MAJOR) OR
	      NOT (STARPU_FIND_VERSION_MINOR STREQUAL STARPU_VERSION_MINOR) )
	    if(STARPU_FIND_REQUIRED AND NOT STARPU_FIND_QUIETLY)
	      message(FATAL_ERROR
		"STARPU version found is ${STARPU_VERSION_STRING} "
		"when required is ${STARPU_FIND_VERSION}")
	    endif()
	  endif()
	endif()
      else()
	if (STARPU_FIND_VERSION_MAJOR AND STARPU_FIND_VERSION_MINOR)
	  # if the version found is older than the required then error
	  if( (STARPU_FIND_VERSION_MAJOR STRGREATER STARPU_VERSION_MAJOR) OR
	      (STARPU_FIND_VERSION_MINOR STRGREATER STARPU_VERSION_MINOR) )
	    if(STARPU_FIND_REQUIRED AND NOT STARPU_FIND_QUIETLY)
	      message(FATAL_ERROR
		"STARPU version found is ${STARPU_VERSION_STRING} "
		"when required is ${STARPU_FIND_VERSION} or newer")
	    endif()
	  endif()
	endif()
      endif()
    else()
      if(STARPU_FIND_REQUIRED AND NOT STARPU_FIND_QUIETLY)
	message(FATAL_ERROR
	  "STARPU version has not been found using starpu_config.h"
	  " located in ${STARPU_starpu_config.h_INCLUDE_DIRS}")
      endif()
    endif()
  else()
    if(STARPU_FIND_REQUIRED AND NOT STARPU_FIND_QUIETLY)
      message(FATAL_ERROR
	"starpu_config.h has not been found while required to get StarPU version")
    endif()
  endif()


  # Try to find the starpu headers in the given paths
  # -------------------------------------------------

  # create list of headers to find
  list(APPEND STARPU_hdrs_to_find "starpu.h;starpu_profiling.h")
  if(STARPU_LOOK_FOR_MPI AND MPI_FOUND)
    list(APPEND STARPU_hdrs_to_find "starpu_mpi.h")
  endif()
  if(STARPU_LOOK_FOR_CUDA AND CUDA_FOUND)
    list(APPEND STARPU_hdrs_to_find "starpu_cuda.h;starpu_scheduler.h")
  endif()

  # call cmake macro to find the header path
  if(STARPU_INCDIR)
    foreach(starpu_hdr ${STARPU_hdrs_to_find})
      set(STARPU_${starpu_hdr}_INCLUDE_DIRS "STARPU_${starpu_hdr}_INCLUDE_DIRS-NOTFOUND")
      find_path(STARPU_${starpu_hdr}_INCLUDE_DIRS
	NAMES ${starpu_hdr}
	HINTS ${STARPU_INCDIR})
    endforeach()
  else()
    if(STARPU_DIR)
      set(STARPU_${starpu_hdr}_INCLUDE_DIRS "STARPU_${starpu_hdr}_INCLUDE_DIRS-NOTFOUND")
      foreach(starpu_hdr ${STARPU_hdrs_to_find})
	find_path(STARPU_${starpu_hdr}_INCLUDE_DIRS
	  NAMES ${starpu_hdr}
	  HINTS ${STARPU_DIR}
	  PATH_SUFFIXES "include"
	  "include/starpu/1.0"
	  "include/starpu/1.1"
	  "include/starpu/1.2"
	  "include/starpu/1.3")
      endforeach()
    else()
      foreach(starpu_hdr ${STARPU_hdrs_to_find})
	set(STARPU_${starpu_hdr}_INCLUDE_DIRS "STARPU_${starpu_hdr}_INCLUDE_DIRS-NOTFOUND")
	find_path(STARPU_${starpu_hdr}_INCLUDE_DIRS
	  NAMES ${starpu_hdr}
	  HINTS ${_inc_env}
	  PATH_SUFFIXES
	  "starpu/1.0"
	  "starpu/1.1"
	  "starpu/1.2"
	  "starpu/1.3")
      endforeach()
    endif()
  endif()

  # If found, add path to cmake variable
  # ------------------------------------
  set(STARPU_INCLUDE_DIRS "")
  foreach(starpu_hdr ${STARPU_hdrs_to_find})

    if (STARPU_${starpu_hdr}_INCLUDE_DIRS)
      # set cmake variables using the pkg-config naming convention
      list(APPEND STARPU_INCLUDE_DIRS "${STARPU_${starpu_hdr}_INCLUDE_DIRS}" )
    else ()
      if(NOT STARPU_FIND_QUIETLY)
	message(STATUS "Looking for starpu -- ${starpu_hdr} not found")
      endif()
      if(starpu_hdr STREQUAL "starpu_mpi.h")
	if(NOT STARPU_FIND_REQUIRED_MPI)
	  if (NOT STARPU_FIND_QUIETLY)
	    message(STATUS "Looking for starpu -- ${starpu_hdr} not required")
	  endif()
	else()
	  list(APPEND STARPU_INCLUDE_DIRS "${STARPU_${starpu_hdr}_INCLUDE_DIRS}" )
	endif()
      elseif( (starpu_hdr STREQUAL "starpu_cuda.h") OR (starpu_hdr STREQUAL "starpu_scheduler.h") )
	if(NOT STARPU_FIND_REQUIRED_CUDA)
	  if (NOT STARPU_FIND_QUIETLY)
	    message(STATUS "Looking for starpu -- ${starpu_hdr} not required")
	  endif()
	else()
	  list(APPEND STARPU_INCLUDE_DIRS "${STARPU_${starpu_hdr}_INCLUDE_DIRS}" )
	endif()
      endif()
    endif ()
    mark_as_advanced(STARPU_${starpu_hdr}_INCLUDE_DIRS)

  endforeach(starpu_hdr ${STARPU_hdrs_to_find})

  if (STARPU_INCLUDE_DIRS)
    list(REMOVE_DUPLICATES STARPU_INCLUDE_DIRS)
  endif ()

  if (STARPU_starpu_config.h_INCLUDE_DIRS)
    # Looking for lib
    # ---------------

    set(STARPU_SHM_LIBRARIES "")
    set(STARPU_MPI_LIBRARIES "")
    set(STARPU_LIBRARY_DIRS "")

    # Add system library paths to search lib
    # --------------------------------------
    unset(_lib_env)
    set(ENV_STARPU_LIBDIR "$ENV{STARPU_LIBDIR}")
    if(ENV_STARPU_LIBDIR)
      list(APPEND _lib_env "${ENV_STARPU_LIBDIR}")
    elseif(ENV_STARPU_DIR)
      list(APPEND _lib_env "${ENV_STARPU_DIR}")
      list(APPEND _lib_env "${ENV_STARPU_DIR}/lib")
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

    # Try to find the starpu libs in the given paths
    # ----------------------------------------------

    # create list of libs to find
    set(STARPU_libs_to_find     "starpu-${STARPU_VERSION_STRING}")
    set(STARPU_SHM_libs_to_find "starpu-${STARPU_VERSION_STRING}")
    if (STARPU_LOOK_FOR_MPI OR MPI_FOUND)
      list(INSERT STARPU_libs_to_find 0 "starpumpi-${STARPU_VERSION_STRING}")
      set(STARPU_MPI_libs_to_find "${STARPU_libs_to_find}")
    endif()

    # call cmake macro to find the lib path
    if(STARPU_LIBDIR)
      foreach(starpu_lib ${STARPU_libs_to_find})
	set(STARPU_${starpu_lib}_LIBRARY "STARPU_${starpu_lib}_LIBRARY-NOTFOUND")
	find_library(STARPU_${starpu_lib}_LIBRARY
	  NAMES ${starpu_lib}
	  HINTS ${STARPU_LIBDIR})
      endforeach()
    else()
      if(STARPU_DIR)
	foreach(starpu_lib ${STARPU_libs_to_find})
	  set(STARPU_${starpu_lib}_LIBRARY "STARPU_${starpu_lib}_LIBRARY-NOTFOUND")
	  find_library(STARPU_${starpu_lib}_LIBRARY
	    NAMES ${starpu_lib}
	    HINTS ${STARPU_DIR}
	    PATH_SUFFIXES lib lib32 lib64)
	endforeach()
      else()
	foreach(starpu_lib ${STARPU_libs_to_find})
	  set(STARPU_${starpu_lib}_LIBRARY "STARPU_${starpu_lib}_LIBRARY-NOTFOUND")
	  find_library(STARPU_${starpu_lib}_LIBRARY
	    NAMES ${starpu_lib}
	    HINTS ${_lib_env})
	endforeach()
      endif()
    endif()

    # If found, add path to cmake variable
    # ------------------------------------
    foreach(starpu_lib ${STARPU_libs_to_find})

      if (STARPU_${starpu_lib}_LIBRARY)

	get_filename_component(${starpu_lib}_lib_path ${STARPU_${starpu_lib}_LIBRARY} PATH)
	# set cmake variables (respects naming convention)

	foreach(starpu_shm_lib ${STARPU_SHM_libs_to_find})
	  if(starpu_shm_lib STREQUAL starpu_lib)
	    list(APPEND STARPU_SHM_LIBRARIES "${STARPU_${starpu_lib}_LIBRARY}")
	  endif()
	endforeach()
	if (STARPU_LOOK_FOR_MPI AND MPI_FOUND)
	  foreach(starpu_mpi_lib ${STARPU_MPI_libs_to_find})
	    if(starpu_mpi_lib STREQUAL starpu_lib)
	      list(APPEND STARPU_MPI_LIBRARIES "${STARPU_${starpu_lib}_LIBRARY}")
	    endif()
	  endforeach()
	endif ()
	list(APPEND STARPU_LIBRARY_DIRS "${${starpu_lib}_lib_path}")

      else (STARPU_${starpu_lib}_LIBRARY)

	if(NOT STARPU_FIND_QUIETLY)
	  message(STATUS "Looking for starpu -- lib ${starpu_lib} not found")
	endif()
	if(starpu_lib STREQUAL "starpumpi-${STARPU_VERSION_STRING}" AND
	    NOT STARPU_FIND_REQUIRED_MPI)
	  # if MPI optional, not a problem: no NOTFOUND in list of MPI LIBRARIES
	  if(NOT STARPU_FIND_QUIETLY)
	    message(STATUS "Looking for starpu -- lib ${starpu_lib} not required")
	  endif()
	else()
	  # for any other lib, add NOTFOUND in the proper list of LIBRARIES
	  foreach(starpu_shm_lib ${STARPU_SHM_libs_to_find})
	    if(starpu_shm_lib STREQUAL starpu_lib)
	      set(STARPU_SHM_LIBRARIES "${STARPU_${starpu_lib}_LIBRARY}")
	    endif()
	  endforeach()
	  if (STARPU_LOOK_FOR_MPI AND MPI_FOUND)
	    foreach(starpu_mpi_lib ${STARPU_MPI_libs_to_find})
	      if(starpu_mpi_lib STREQUAL starpu_lib)
		list(APPEND STARPU_MPI_LIBRARIES "${STARPU_${starpu_mpi_lib}_LIBRARY}")
	      endif()
	    endforeach()
	  endif ()
	endif()

      endif (STARPU_${starpu_lib}_LIBRARY)

      mark_as_advanced(STARPU_${starpu_lib}_LIBRARY)

    endforeach(starpu_lib ${STARPU_libs_to_find})

    list(REMOVE_DUPLICATES STARPU_SHM_LIBRARIES)
    list(REMOVE_DUPLICATES STARPU_MPI_LIBRARIES)
    if (STARPU_LIBRARY_DIRS)
      list(REMOVE_DUPLICATES STARPU_LIBRARY_DIRS)
    endif ()

    if (STARPU_SHM_LIBRARIES AND STARPU_LIBRARY_DIRS AND STARPU_INCLUDE_DIRS)
      set(STARPU_SHM_FOUND TRUE)
      if(STARPU_MPI_LIBRARIES)
	set(STARPU_MPI_FOUND TRUE)
      endif()
    endif()

  else(STARPU_starpu_config.h_INCLUDE_DIRS)
    if(NOT STARPU_FIND_QUIETLY)
      message(STATUS "The version of StarPU is not known so that we do not search libraries")
    endif()
  endif(STARPU_starpu_config.h_INCLUDE_DIRS)

  foreach(lib ${STARPU_SHM_LIBRARIES})
    if (NOT lib)
      set(STARPU_SHM_LIBRARIES "STARPU_SHM_LIBRARIES-NOTFOUND")
    endif()
  endforeach()
  foreach(lib ${STARPU_MPI_LIBRARIES})
    if (NOT lib)
      set(STARPU_MPI_LIBRARIES "STARPU_MPI_LIBRARIES-NOTFOUND")
    endif()
  endforeach()
  if(STARPU_LOOK_FOR_MPI)
    set(STARPU_LIBRARIES "${STARPU_MPI_LIBRARIES}")
  else()
    set(STARPU_LIBRARIES "${STARPU_SHM_LIBRARIES}")
  endif()

  # check a function to validate the find
  if(STARPU_LIBRARIES)

    set(REQUIRED_FLAGS)
    set(REQUIRED_LDFLAGS)
    set(REQUIRED_INCDIRS)
    set(REQUIRED_LIBDIRS)
    set(REQUIRED_LIBS)

    # STARPU
    if (STARPU_INCLUDE_DIRS)
      set(REQUIRED_INCDIRS "${STARPU_INCLUDE_DIRS}")
    endif()
    set(CMAKE_REQUIRED_FLAGS)
    foreach(libdir ${STARPU_LIBRARY_DIRS})
      if (libdir)
	list(APPEND REQUIRED_LIBDIRS "${libdir}")
      endif()
    endforeach()
    set(REQUIRED_LIBS "${STARPU_LIBRARIES}")
    # HWLOC
    if (HWLOC_FOUND AND STARPU_LOOK_FOR_HWLOC)
      if (HWLOC_INCLUDE_DIRS)
	list(APPEND REQUIRED_INCDIRS "${HWLOC_INCLUDE_DIRS}")
      endif()
      if (HWLOC_LIBRARY_DIRS)
	list(APPEND REQUIRED_LIBDIRS "${HWLOC_LIBRARY_DIRS}")
      endif()
      foreach(lib ${HWLOC_LIBRARIES})
	if (EXISTS ${lib} OR ${lib} MATCHES "^-")
	  list(APPEND REQUIRED_LIBS "${lib}")
	else()
	  list(APPEND REQUIRED_LIBS "-l${lib}")
	endif()
      endforeach()
    endif()
    # MPI
    if (MPI_FOUND AND STARPU_LOOK_FOR_MPI)
      if (MPI_C_INCLUDE_PATH)
	list(APPEND REQUIRED_INCDIRS "${MPI_C_INCLUDE_PATH}")
      endif()
      if (MPI_C_LINK_FLAGS)
	if (${MPI_C_LINK_FLAGS} MATCHES "  -")
	  string(REGEX REPLACE " -" "-" MPI_C_LINK_FLAGS ${MPI_C_LINK_FLAGS})
	endif()
	list(APPEND REQUIRED_LDFLAGS "${MPI_C_LINK_FLAGS}")
      endif()
      list(APPEND REQUIRED_LIBS "${MPI_C_LIBRARIES}")
    endif()
    # MAGMA
    if (MAGMA_FOUND AND STARPU_LOOK_FOR_MAGMA)
      if (MAGMA_INCLUDE_DIRS_DEP)
	list(APPEND REQUIRED_INCDIRS "${MAGMA_INCLUDE_DIRS_DEP}")
      elseif(MAGMA_INCLUDE_DIRS)
	list(APPEND REQUIRED_INCDIRS "${MAGMA_INCLUDE_DIRS}")
      endif()
      if (MAGMA_LIBRARY_DIRS_DEP)
	list(APPEND REQUIRED_LIBDIRS "${MAGMA_LIBRARY_DIRS_DEP}")
      elseif(MAGMA_LIBRARY_DIRS)
	list(APPEND REQUIRED_LIBDIRS "${MAGMA_LIBRARY_DIRS}")
      endif()
      if (MAGMA_LIBRARIES_DEP)
	list(APPEND REQUIRED_LIBS "${MAGMA_LIBRARIES_DEP}")
      elseif(MAGMA_LIBRARIES)
	foreach(lib ${MAGMA_LIBRARIES})
	  if (EXISTS ${lib} OR ${lib} MATCHES "^-")
	    list(APPEND REQUIRED_LIBS "${lib}")
	  else()
	    list(APPEND REQUIRED_LIBS "-l${lib}")
	  endif()
	endforeach()
      endif()
      if (MAGMA_LINKER_FLAGS)
	list(APPEND REQUIRED_LDFLAGS "${MAGMA_LINKER_FLAGS}")
      endif()
    endif()
    # CUDA
    if (CUDA_FOUND AND STARPU_LOOK_FOR_CUDA)
      if (CUDA_INCLUDE_DIRS)
	list(APPEND REQUIRED_INCDIRS "${CUDA_INCLUDE_DIRS}")
      endif()
      if (CUDA_LIBRARY_DIRS)
	list(APPEND REQUIRED_LIBDIRS "${CUDA_LIBRARY_DIRS}")
      endif()
      list(APPEND REQUIRED_LIBS "${CUDA_CUBLAS_LIBRARIES};${CUDA_LIBRARIES}")
    endif()
    # FXT
    if (FXT_FOUND AND STARPU_LOOK_FOR_FXT)
      if (FXT_INCLUDE_DIRS)
	list(APPEND REQUIRED_INCDIRS "${FXT_INCLUDE_DIRS}")
      endif()
      if (FXT_LIBRARY_DIRS)
	list(APPEND REQUIRED_LIBDIRS "${FXT_LIBRARY_DIRS}")
      endif()
      foreach(lib ${FXT_LIBRARIES})
	if (EXISTS ${lib} OR ${lib} MATCHES "^-")
	  list(APPEND REQUIRED_LIBS "${lib}")
	else()
	  list(APPEND REQUIRED_LIBS "-l${lib}")
	endif()
      endforeach()
    endif()
    # SIMGRID
    if (SIMGRID_FOUND AND STARPU_LOOK_FOR_SIMGRID)
      if (SIMGRID_INCLUDE_DIRS)
	list(APPEND REQUIRED_INCDIRS "${SIMGRID_INCLUDE_DIRS}")
      endif()
      if (SIMGRID_LIBRARY_DIRS)
	list(APPEND REQUIRED_LIBDIRS "${SIMGRID_LIBRARY_DIRS}")
      endif()
      foreach(lib ${SIMGRID_LIBRARIES})
	if (EXISTS ${lib} OR ${lib} MATCHES "^-")
	  list(APPEND REQUIRED_LIBS "${lib}")
	else()
	  list(APPEND REQUIRED_LIBS "-l${lib}")
	endif()
      endforeach()
      list(APPEND REQUIRED_FLAGS "-include starpu_simgrid_wrap.h")
    endif()
    # BLAS
    if (BLAS_FOUND AND STARPU_LOOK_FOR_BLAS)
      if (BLAS_INCLUDE_DIRS)
	list(APPEND REQUIRED_INCDIRS "${BLAS_INCLUDE_DIRS}")
      endif()
      if (BLAS_LIBRARY_DIRS)
	list(APPEND REQUIRED_LIBDIRS "${BLAS_LIBRARY_DIRS}")
      endif()
      list(APPEND REQUIRED_LIBS "${BLAS_LIBRARIES}")
      if (BLAS_LINKER_FLAGS)
	list(APPEND REQUIRED_LDFLAGS "${BLAS_LINKER_FLAGS}")
      endif()
    endif()
    # Fortran
    if (CMAKE_C_COMPILER_ID MATCHES "GNU")
      find_library(
	FORTRAN_gfortran_LIBRARY
	NAMES gfortran
	HINTS ${_lib_env}
	)
      mark_as_advanced(FORTRAN_gfortran_LIBRARY)
      if (FORTRAN_gfortran_LIBRARY AND CMAKE_C_COMPILER_ID STREQUAL "GNU")
	list(APPEND REQUIRED_LIBS "${FORTRAN_gfortran_LIBRARY}")
      endif()
    elseif (CMAKE_C_COMPILER_ID MATCHES "Intel")
      find_library(
	FORTRAN_ifcore_LIBRARY
	NAMES ifcore
	HINTS ${_lib_env}
	)
      mark_as_advanced(FORTRAN_ifcore_LIBRARY)
      if (FORTRAN_ifcore_LIBRARY)
	list(APPEND REQUIRED_LIBS "${FORTRAN_ifcore_LIBRARY}")
      endif()
    endif()
    # EXTRA LIBS such that pthread, m, rt
    list(APPEND REQUIRED_LIBS ${STARPU_EXTRA_LIBRARIES})

    # set required libraries for link
    set(CMAKE_REQUIRED_INCLUDES "${REQUIRED_INCDIRS}")
    set(CMAKE_REQUIRED_LIBRARIES)
    list(APPEND CMAKE_REQUIRED_LIBRARIES "${REQUIRED_LDFLAGS}")
    foreach(lib_dir ${REQUIRED_LIBDIRS})
      list(APPEND CMAKE_REQUIRED_LIBRARIES "-L${lib_dir}")
    endforeach()
    list(APPEND CMAKE_REQUIRED_LIBRARIES "${REQUIRED_LIBS}")
    list(APPEND CMAKE_REQUIRED_FLAGS "${REQUIRED_FLAGS}")
    string(REGEX REPLACE "^ -" "-" CMAKE_REQUIRED_LIBRARIES "${CMAKE_REQUIRED_LIBRARIES}")

    # test link
    unset(STARPU_WORKS CACHE)
    include(CheckFunctionExists)
    check_function_exists(starpu_init STARPU_WORKS)
    mark_as_advanced(STARPU_WORKS)

    if(STARPU_WORKS)
      # save link with dependencies
      if (REQUIRED_FLAGS)
	set(STARPU_LIBRARIES_DEP "${REQUIRED_FLAGS};${REQUIRED_LIBS}")
      else()
	set(STARPU_LIBRARIES_DEP "${REQUIRED_LIBS}")
      endif()
      set(STARPU_LIBRARY_DIRS_DEP "${REQUIRED_LIBDIRS}")
      set(STARPU_INCLUDE_DIRS_DEP "${REQUIRED_INCDIRS}")
      set(STARPU_LINKER_FLAGS "${REQUIRED_LDFLAGS}")
      list(REMOVE_DUPLICATES STARPU_LIBRARY_DIRS_DEP)
      list(REMOVE_DUPLICATES STARPU_INCLUDE_DIRS_DEP)
      list(REMOVE_DUPLICATES STARPU_LINKER_FLAGS)
    else()
      if(NOT STARPU_FIND_QUIETLY)
	message(STATUS "Looking for starpu : test of starpu_init fails")
	message(STATUS "CMAKE_REQUIRED_LIBRARIES: ${CMAKE_REQUIRED_LIBRARIES}")
	message(STATUS "CMAKE_REQUIRED_INCLUDES: ${CMAKE_REQUIRED_INCLUDES}")
	message(STATUS "Check in CMakeFiles/CMakeError.log to figure out why it fails")
	message(STATUS "Maybe STARPU is linked with specific libraries. "
	  "Have you tried with COMPONENTS (HWLOC, CUDA, MPI, BLAS, MAGMA, FXT, SIMGRID)? "
	  "See the explanation in FindSTARPU.cmake.")
      endif()
    endif()
    set(CMAKE_REQUIRED_INCLUDES)
    set(CMAKE_REQUIRED_FLAGS)
    set(CMAKE_REQUIRED_LIBRARIES)
  endif(STARPU_LIBRARIES)

endif( (NOT PKG_CONFIG_EXECUTABLE) OR (PKG_CONFIG_EXECUTABLE AND NOT STARPU_FOUND) OR (STARPU_GIVEN_BY_USER) )

if (STARPU_LIBRARIES)
  if (STARPU_LIBRARY_DIRS)
    foreach(dir ${STARPU_LIBRARY_DIRS})
      if ("${dir}" MATCHES "starpu")
	set(first_lib_path "${dir}")
      endif()
    endforeach()
  else()
    list(GET STARPU_LIBRARIES 0 first_lib)
    get_filename_component(first_lib_path "${first_lib}" PATH)
  endif()
  if (${first_lib_path} MATCHES "/lib(32|64)?$")
    string(REGEX REPLACE "/lib(32|64)?$" "" not_cached_dir "${first_lib_path}")
    set(STARPU_DIR_FOUND "${not_cached_dir}" CACHE PATH "Installation directory of STARPU library" FORCE)
  else()
    set(STARPU_DIR_FOUND "${first_lib_path}" CACHE PATH "Installation directory of STARPU library" FORCE)
  endif()
endif()
mark_as_advanced(STARPU_DIR)
mark_as_advanced(STARPU_DIR_FOUND)

# check that STARPU has been found
# --------------------------------
include(FindPackageHandleStandardArgs)
if(NOT STARPU_FIND_QUIETLY)
  if(STARPU_SHM_FOUND)
    message(STATUS "StarPU library has been found.")
    if(STARPU_MPI_LIBRARIES)
      message(STATUS "The mpi version of StarPU has been found so that we manage"
	" two lists of libs, one sequential and one parallel (see"
	" STARPU_SHM_LIBRARIES and STARPU_MPI_LIBRARIES).")
    endif()
    message(STATUS "StarPU shared memory libraries stored in STARPU_SHM_LIBRARIES")
  endif()
endif()
if (PKG_CONFIG_EXECUTABLE AND STARPU_SHM_FOUND AND NOT STARPU_GIVEN_BY_USER)
  find_package_handle_standard_args(STARPU DEFAULT_MSG
    STARPU_SHM_LIBRARIES)
else()
  find_package_handle_standard_args(STARPU DEFAULT_MSG
    STARPU_SHM_LIBRARIES
    STARPU_WORKS)
endif()
if(STARPU_LOOK_FOR_MPI)
  if(STARPU_MPI_LIBRARIES AND NOT STARPU_FIND_QUIETLY)
    message(STATUS "StarPU mpi libraries stored in STARPU_MPI_LIBRARIES")
  endif()
  if (PKG_CONFIG_EXECUTABLE AND STARPU_MPI_FOUND AND NOT STARPU_GIVEN_BY_USER)
    find_package_handle_standard_args(STARPU DEFAULT_MSG
      STARPU_MPI_LIBRARIES)
  else()
    find_package_handle_standard_args(STARPU DEFAULT_MSG
      STARPU_MPI_LIBRARIES
      STARPU_WORKS)
  endif()
endif()
