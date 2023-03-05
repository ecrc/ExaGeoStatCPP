###
#
# @copyright (c) 2009-2014 The University of Tennessee and The University
#                          of Tennessee Research Foundation.
#                          All rights reserved.
# @copyright (c) 2012-2016 Inria. All rights reserved.
# @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
#
# @copyright (c) 2022 King Abdullah University of Science and Technology (KAUST).
#                     All rights reserved.
#
###
#
# - Find AL4SAN include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(AL4SAN
#               [REQUIRED]             # Fail with error if al4san is not found
#               [COMPONENTS <comp1> <comp2> ...] # dependencies
#              )
#
#  AL4SAN depends on the following libraries:
#   - Threads, m, rt
#   - HWLOC
#   - CBLAS
#   - LAPACKE
#   - TMG
#   - At least one runtime, default is StarPU
#     (For QUARK, use COMPONENTS QUARK)
#
#  COMPONENTS are optional libraries AL4SAN could be linked with,
#  Use it to drive detection of a specific compilation chain
#  COMPONENTS can be some of the following:
#   - STARPU (default): to activate detection of AL4san linked with StarPU
#   - QUARK (STARPU will be deactivated): to activate detection of AL4san linked with QUARK
#   - CUDA (comes with cuBLAS): to activate detection of AL4san linked with CUDA
#   - MAGMA: to activate detection of AL4san linked with MAGMA
#   - MPI: to activate detection of AL4san linked with MPI
#   - FXT: to activate detection of AL4san linked with StarPU+FXT
#
# This module finds headers and al4san library.
# Results are reported in variables:
#  AL4SAN_FOUND            - True if headers and requested libraries were found
#  AL4SAN_C_FLAGS          - list of required compilation flags (excluding -I)
#  AL4SAN_LINKER_FLAGS     - list of required linker flags (excluding -l and -L)
#  AL4SAN_INCLUDE_DIRS     - al4san include directories
#  AL4SAN_LIBRARY_DIRS     - Link directories for al4san libraries
#  AL4SAN_INCLUDE_DIRS_DEP - al4san + dependencies include directories
#  AL4SAN_LIBRARY_DIRS_DEP - al4san + dependencies link directories
#  AL4SAN_LIBRARIES_DEP    - al4san libraries + dependencies
# The user can give specific paths where to find the libraries adding cmake
# options at configure (ex: cmake path/to/project -DAL4SAN_DIR=path/to/al4san):
#  AL4SAN_DIR              - Where to find the base directory of al4san
#  AL4SAN_INCDIR           - Where to find the header files
#  AL4SAN_LIBDIR           - Where to find the library files
# The module can also look for the following environment variables if paths
# are not given as cmake variable: AL4SAN_DIR, AL4SAN_INCDIR, AL4SAN_LIBDIR

#=============================================================================
# Copyright 2012-2013 Inria
# Copyright 2012-2013 Emmanuel Agullo
# Copyright 2012-2013 Mathieu Faverge
# Copyright 2012      Cedric Castagnede
# Copyright 2013-2016 Florent Pruvost
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file AL4SAN-Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of AL4san, substitute the full
#  License text for the above reference.)


if (NOT AL4SAN_FOUND)
  set(AL4SAN_DIR "" CACHE PATH "Installation directory of AL4SAN library")
  if (NOT AL4SAN_FIND_QUIETLY)
    message(STATUS "A cache variable, namely AL4SAN_DIR, has been set to specify the install directory of AL4SAN")
  endif()
endif()

# Try to find AL4SAN dependencies if specified as COMPONENTS during the call
set(AL4SAN_LOOK_FOR_STARPU ON)
set(AL4SAN_LOOK_FOR_QUARK OFF)
set(AL4SAN_LOOK_FOR_OPENMP OFF)
set(AL4SAN_LOOK_FOR_PARSEC OFF)
set(AL4SAN_LOOK_FOR_CUDA OFF)
set(AL4SAN_LOOK_FOR_MAGMA OFF)
set(AL4SAN_LOOK_FOR_MPI OFF)
set(AL4SAN_LOOK_FOR_FXT OFF)

if( AL4SAN_FIND_COMPONENTS )
  foreach( component ${AL4SAN_FIND_COMPONENTS} )
    if (${component} STREQUAL "STARPU")
      # means we look for AL4san with StarPU
      set(AL4SAN_LOOK_FOR_STARPU ON)
      set(AL4SAN_LOOK_FOR_QUARK OFF)
      set(AL4SAN_LOOK_FOR_OPENMP OFF)
      set(AL4SAN_LOOK_FOR_PARSEC OFF)
    endif()
    if (${component} STREQUAL "QUARK")
      # means we look for AL4san with QUARK
      set(AL4SAN_LOOK_FOR_QUARK ON)
      set(AL4SAN_LOOK_FOR_STARPU OFF)
      set(AL4SAN_LOOK_FOR_PARSEC OFF)
      set(AL4SAN_LOOK_FOR_OPENMP OFF)
    endif()
     if (${component} STREQUAL "PARSEC")
      # means we look for AL4san with PARSEC
       set(AL4SAN_LOOK_FOR_PARSEC ON)
       set(AL4SAN_LOOK_FOR_STARPU OFF)
       set(AL4SAN_LOOK_FOR_QUARK OFF)
      set(AL4SAN_LOOK_FOR_OPENMP OFF)
     endif()
     if (${component} STREQUAL "OPENMP")
      # means we look for AL4san with PARSEC
       set(AL4SAN_LOOK_FOR_OPENMP ON)
       set(AL4SAN_LOOK_FOR_STARPU OFF)
       set(AL4SAN_LOOK_FOR_QUARK OFF)
      set(AL4SAN_LOOK_FOR_PARSEC OFF)
     endif()
    if (${component} STREQUAL "CUDA")
      # means we look for AL4san with CUDA
      set(AL4SAN_LOOK_FOR_CUDA ON)
    endif()
    if (${component} STREQUAL "MAGMA")
      # means we look for AL4san with MAGMA
      set(AL4SAN_LOOK_FOR_MAGMA ON)
    endif()
    if (${component} STREQUAL "MPI")
      # means we look for AL4san with MPI
      set(AL4SAN_LOOK_FOR_MPI ON)
    endif()
    if (${component} STREQUAL "FXT")
      # means we look for AL4san with FXT
      set(AL4SAN_LOOK_FOR_FXT ON)
    endif()
  endforeach()
endif()

set(ENV_AL4SAN_DIR "$ENV{AL4SAN_DIR}")
set(ENV_AL4SAN_INCDIR "$ENV{AL4SAN_INCDIR}")
set(ENV_AL4SAN_LIBDIR "$ENV{AL4SAN_LIBDIR}")
set(AL4SAN_GIVEN_BY_USER "FALSE")
if ( AL4SAN_DIR OR ( AL4SAN_INCDIR AND AL4SAN_LIBDIR) OR ENV_AL4SAN_DIR OR (ENV_AL4SAN_INCDIR AND ENV_AL4SAN_LIBDIR) )
  set(AL4SAN_GIVEN_BY_USER "TRUE")
endif()

# Optionally use pkg-config to detect include/library dirs (if pkg-config is available)
# -------------------------------------------------------------------------------------
include(FindPkgConfig)
find_package(PkgConfig QUIET)
if(PKG_CONFIG_EXECUTABLE AND NOT AL4SAN_GIVEN_BY_USER)

  pkg_search_module(AL4SAN al4san)
  if (NOT AL4SAN_FIND_QUIETLY)
    if (AL4SAN_FOUND AND AL4SAN_LIBRARIES)
      message(STATUS "Looking for AL4SAN - found using PkgConfig")
      #if(NOT AL4SAN_INCLUDE_DIRS)
      #    message("${Magenta}AL4SAN_INCLUDE_DIRS is empty using PkgConfig."
      #        "Perhaps the path to al4san headers is already present in your"
      #        "C(PLUS)_INCLUDE_PATH environment variable.${ColourReset}")
      #endif()
    else()
      message(STATUS "${Magenta}Looking for AL4SAN - not found using PkgConfig."
	"\n   Perhaps you should add the directory containing al4san.pc"
	"\n   to the PKG_CONFIG_PATH environment variable.${ColourReset}")
    endif()
  endif()

  if (AL4SAN_FIND_VERSION_EXACT)
    if( NOT (AL4SAN_FIND_VERSION_MAJOR STREQUAL AL4SAN_VERSION_MAJOR) OR
	NOT (AL4SAN_FIND_VERSION_MINOR STREQUAL AL4SAN_VERSION_MINOR) )
      if(NOT AL4SAN_FIND_QUIETLY)
	message(FATAL_ERROR
	  "AL4SAN version found is ${AL4SAN_VERSION_STRING}"
	  "when required is ${AL4SAN_FIND_VERSION}")
      endif()
    endif()
  else()
    # if the version found is older than the required then error
    if( (AL4SAN_FIND_VERSION_MAJOR STRGREATER AL4SAN_VERSION_MAJOR) OR
	(AL4SAN_FIND_VERSION_MINOR STRGREATER AL4SAN_VERSION_MINOR) )
      if(NOT AL4SAN_FIND_QUIETLY)
	message(FATAL_ERROR
	  "AL4SAN version found is ${AL4SAN_VERSION_STRING}"
	  "when required is ${AL4SAN_FIND_VERSION} or newer")
      endif()
    endif()
  endif()

  set(AL4SAN_C_FLAGS "${AL4SAN_CFLAGS_OTHER}")
  set(AL4SAN_INCLUDE_DIRS_DEP "${AL4SAN_STATIC_INCLUDE_DIRS}")
  set(AL4SAN_LIBRARY_DIRS_DEP "${AL4SAN_STATIC_LIBRARY_DIRS}")
  set(AL4SAN_LIBRARIES_DEP "${AL4SAN_STATIC_LIBRARIES}")

endif(PKG_CONFIG_EXECUTABLE AND NOT AL4SAN_GIVEN_BY_USER)

if( (NOT PKG_CONFIG_EXECUTABLE) OR (PKG_CONFIG_EXECUTABLE AND NOT AL4SAN_FOUND) OR (AL4SAN_GIVEN_BY_USER) )

  if (NOT AL4SAN_FIND_QUIETLY)
    message(STATUS "Looking for AL4SAN - PkgConfig not used")
  endif()

  # Dependencies detection
  # ----------------------

  if (NOT AL4SAN_FIND_QUIETLY)
    message(STATUS "Looking for AL4SAN - Try to detect pthread")
  endif()
  if (AL4SAN_FIND_REQUIRED)
    find_package(Threads REQUIRED)
  else()
    find_package(Threads)
  endif()
  set(AL4SAN_EXTRA_LIBRARIES "")
  if( THREADS_FOUND )
    list(APPEND AL4SAN_EXTRA_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})
  endif ()

  # Add math library to the list of extra
  # it normally exists on all common systems provided with a C compiler
  if (NOT AL4SAN_FIND_QUIETLY)
    message(STATUS "Looking for AL4SAN - Try to detect libm")
  endif()
  set(AL4SAN_M_LIBRARIES "")
  if(UNIX OR WIN32)
    find_library(
      AL4SAN_M_m_LIBRARY
      NAMES m
      )
    mark_as_advanced(AL4SAN_M_m_LIBRARY)
    if (AL4SAN_M_m_LIBRARY)
      list(APPEND AL4SAN_M_LIBRARIES "${AL4SAN_M_m_LIBRARY}")
      list(APPEND AL4SAN_EXTRA_LIBRARIES "${AL4SAN_M_m_LIBRARY}")
    else()
      if (AL4SAN_FIND_REQUIRED)
	message(FATAL_ERROR "Could NOT find libm on your system."
	  "Are you sure to a have a C compiler installed?")
      endif()
    endif()
  endif()

  # Try to find librt (libposix4 - POSIX.1b Realtime Extensions library)
  # on Unix systems except Apple ones because it does not exist on it
  if (NOT AL4SAN_FIND_QUIETLY)
    message(STATUS "Looking for AL4SAN - Try to detect librt")
  endif()
  set(AL4SAN_RT_LIBRARIES "")
  if(UNIX AND NOT APPLE)
    find_library(
      AL4SAN_RT_rt_LIBRARY
      NAMES rt
      )
    mark_as_advanced(AL4SAN_RT_rt_LIBRARY)
    if (AL4SAN_RT_rt_LIBRARY)
      list(APPEND AL4SAN_RT_LIBRARIES "${AL4SAN_RT_rt_LIBRARY}")
      list(APPEND AL4SAN_EXTRA_LIBRARIES "${AL4SAN_RT_rt_LIBRARY}")
    else()
      if (AL4SAN_FIND_REQUIRED)
	message(FATAL_ERROR "Could NOT find librt on your system")
      endif()
    endif()
  endif()

  # AL4SAN depends on CBLAS
  #---------------------------
  if (NOT AL4SAN_FIND_QUIETLY)
    message(STATUS "Looking for AL4SAN - Try to detect CBLAS (depends on BLAS)")
  endif()
  if (AL4SAN_FIND_REQUIRED)
    find_package(CBLAS REQUIRED)
  else()
    find_package(CBLAS)
  endif()

  # AL4SAN depends on LAPACKE
  #-----------------------------

  # standalone version of lapacke seems useless for now
  # let the comment in case we meet some problems of non existing lapacke
  # functions in lapack library such as mkl, acml, ...
  #set(LAPACKE_STANDALONE TRUE)
  if (NOT AL4SAN_FIND_QUIETLY)
    message(STATUS "Looking for AL4SAN - Try to detect LAPACKE (depends on LAPACK)")
  endif()
  if (AL4SAN_FIND_REQUIRED)
    find_package(LAPACKE REQUIRED)
  else()
    find_package(LAPACKE)
  endif()

  # AL4SAN depends on TMG
  #-------------------------
  if (NOT AL4SAN_FIND_QUIETLY)
    message(STATUS "Looking for AL4SAN - Try to detect TMG (depends on LAPACK)")
  endif()
  if (AL4SAN_FIND_REQUIRED)
    find_package(TMG REQUIRED)
  else()
    find_package(TMG)
  endif()

  # AL4SAN may depend on CUDA/CUBLAS
  #------------------------------------
  if (NOT CUDA_FOUND AND AL4SAN_LOOK_FOR_CUDA)
    if (AL4SAN_FIND_REQUIRED AND AL4SAN_FIND_REQUIRED_CUDA)
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

  # AL4SAN may depend on MAGMA gpu kernels
  # call our cmake module to test (in cmake_modules)
  # change this call position if not appropriated
  #-------------------------------------------------
  if( CUDA_FOUND AND AL4SAN_LOOK_FOR_MAGMA )
    set(AL4SAN_MAGMA_VERSION "1.4" CACHE STRING "oldest MAGMA version desired")
    if (AL4SAN_FIND_REQUIRED AND AL4SAN_FIND_REQUIRED_MAGMA)
      find_package(MAGMA ${AL4SAN_MAGMA_VERSION} REQUIRED)
    else()
      find_package(MAGMA ${AL4SAN_MAGMA_VERSION})
    endif()
  endif()

  # AL4SAN depends on MPI
  #-------------------------
  if( NOT MPI_FOUND AND AL4SAN_LOOK_FOR_MPI )

    # allows to use an external mpi compilation by setting compilers with
    # -DMPI_C_COMPILER=path/to/mpicc -DMPI_Fortran_COMPILER=path/to/mpif90
    # at cmake configure
    if(NOT MPI_C_COMPILER)
      set(MPI_C_COMPILER mpicc)
    endif()
    if (AL4SAN_FIND_REQUIRED AND AL4SAN_FIND_REQUIRED_MPI)
      find_package(MPI REQUIRED)
    else()
      find_package(MPI)
    endif()
    if (MPI_FOUND)
      mark_as_advanced(MPI_LIBRARY)
      mark_as_advanced(MPI_EXTRA_LIBRARY)
    endif()

  endif()

  if( NOT STARPU_FOUND AND AL4SAN_LOOK_FOR_STARPU )

    set(AL4SAN_STARPU_VERSION "1.1" CACHE STRING "oldest STARPU version desired")

    # create list of components in order to make a single call to find_package(starpu...)
    # we explicitly need a StarPU version built with hwloc
    set(STARPU_COMPONENT_LIST "HWLOC")

    # StarPU may depend on MPI
    # allows to use an external mpi compilation by setting compilers with
    # -DMPI_C_COMPILER=path/to/mpicc -DMPI_Fortran_COMPILER=path/to/mpif90
    # at cmake configure
    if (AL4SAN_LOOK_FOR_MPI)
      if(NOT MPI_C_COMPILER)
	set(MPI_C_COMPILER mpicc)
      endif()
      list(APPEND STARPU_COMPONENT_LIST "MPI")
    endif()
    if (AL4SAN_LOOK_FOR_CUDA)
      list(APPEND STARPU_COMPONENT_LIST "CUDA")
    endif()
    if (AL4SAN_LOOK_FOR_FXT)
      list(APPEND STARPU_COMPONENT_LIST "FXT")
    endif()
    if (AL4SAN_FIND_REQUIRED AND AL4SAN_FIND_REQUIRED_STARPU)
      find_package(STARPU ${AL4SAN_STARPU_VERSION} REQUIRED
	COMPONENTS ${STARPU_COMPONENT_LIST})
    else()
      find_package(STARPU ${AL4SAN_STARPU_VERSION}
	COMPONENTS ${STARPU_COMPONENT_LIST})
    endif()

  endif()

  if( NOT QUARK_FOUND AND AL4SAN_LOOK_FOR_QUARK )

    # try to find quark runtime
    if (AL4SAN_FIND_REQUIRED AND AL4SAN_FIND_REQUIRED_QUARK)
      find_package(QUARK REQUIRED COMPONENTS HWLOC)
    else()
      find_package(QUARK COMPONENTS HWLOC)
    endif()

  if( NOT PARSEC_FOUND AND AL4SAN_LOOK_FOR_PARSEC)


    # create list of components in order to make a single call to find_package(parsec...)
    # we explicitly need a StarPU version built with hwloc
    set(PARSEC_COMPONENT_LIST "HWLOC")

    # PARSEC may depend on MPI
    # allows to use an external mpi compilation by setting compilers with
    # -DMPI_C_COMPILER=path/to/mpicc -DMPI_Fortran_COMPILER=path/to/mpif90
    # at cmake configure
    if (AL4SAN_LOOK_FOR_MPI)
      if(NOT MPI_C_COMPILER)
        set(MPI_C_COMPILER mpicc)
      endif()
      list(APPEND PARSEC_COMPONENT_LIST "MPI")
    endif()
    if (AL4SAN_LOOK_FOR_CUDA)
      list(APPEND PARSEC_COMPONENT_LIST "CUDA")
    endif()
    if (AL4SAN_FIND_REQUIRED AND AL4SAN_FIND_REQUIRED_PARSEC)
      find_package(PARSEC REQUIRED ${PARSEC_COMPONENT_LIST})
    else()
      find_package(PARSEC COMPONENTS ${PARSEC_COMPONENT_LIST})
    endif()
  endif()

  if( NOT OPENMP_FOUND AND AL4SAN_LOOK_FOR_OPENMP )

    # try to find openmp runtime
    if (AL4SAN_FIND_REQUIRED AND AL4SAN_FIND_REQUIRED_OPENMP)
      find_package(OpenMP)
    else()
      find_package(OpenMP)
    endif()
  endif()

  endif()

  # Looking for include
  # -------------------

  # Add system include paths to search include
  # ------------------------------------------
  unset(_inc_env)
  set(ENV_AL4SAN_DIR "$ENV{AL4SAN_DIR}")
  set(ENV_AL4SAN_INCDIR "$ENV{AL4SAN_INCDIR}")
  if(ENV_AL4SAN_INCDIR)
    list(APPEND _inc_env "${ENV_AL4SAN_INCDIR}")
  elseif(ENV_AL4SAN_DIR)
    list(APPEND _inc_env "${ENV_AL4SAN_DIR}")
    list(APPEND _inc_env "${ENV_AL4SAN_DIR}/include")
    list(APPEND _inc_env "${ENV_AL4SAN_DIR}/include/al4san")
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


  # Try to find the al4san header in the given paths
  # ---------------------------------------------------
  # call cmake macro to find the header path
  if(AL4SAN_INCDIR)
    set(AL4SAN_al4san.h_DIRS "AL4SAN_al4san.h_DIRS-NOTFOUND")
    find_path(AL4SAN_al4san.h_DIRS
      NAMES al4san.h
      HINTS ${AL4SAN_INCDIR})
  else()
    if(AL4SAN_DIR)
      set(AL4SAN_al4san.h_DIRS "AL4SAN_al4san.h_DIRS-NOTFOUND")
      find_path(AL4SAN_al4san.h_DIRS
	NAMES al4san.h
	HINTS ${AL4SAN_DIR}
	PATH_SUFFIXES "include" "include/al4san")
    else()
      set(AL4SAN_al4san.h_DIRS "AL4SAN_al4san.h_DIRS-NOTFOUND")
      find_path(AL4SAN_al4san.h_DIRS
	NAMES al4san.h
	HINTS ${_inc_env}
	PATH_SUFFIXES "al4san")
    endif()
  endif()
  mark_as_advanced(AL4SAN_al4san.h_DIRS)

  # If found, add path to cmake variable
  # ------------------------------------
  if (AL4SAN_al4san.h_DIRS)
    set(AL4SAN_INCLUDE_DIRS "${AL4SAN_al4san.h_DIRS}")
  else ()
    set(AL4SAN_INCLUDE_DIRS "AL4SAN_INCLUDE_DIRS-NOTFOUND")
    if(NOT AL4SAN_FIND_QUIETLY)
      message(STATUS "Looking for al4san -- al4san.h not found")
    endif()
  endif()


  # Looking for lib
  # ---------------

  # Add system library paths to search lib
  # --------------------------------------
  unset(_lib_env)
  set(ENV_AL4SAN_LIBDIR "$ENV{AL4SAN_LIBDIR}")
  if(ENV_AL4SAN_LIBDIR)
    list(APPEND _lib_env "${ENV_AL4SAN_LIBDIR}")
  elseif(ENV_AL4SAN_DIR)
    list(APPEND _lib_env "${ENV_AL4SAN_DIR}")
    list(APPEND _lib_env "${ENV_AL4SAN_DIR}/lib")
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

  # Try to find the al4san lib in the given paths
  # ------------------------------------------------

  # create list of libs to find
  set(AL4SAN_libs_to_find "al4san")
  if (STARPU_FOUND)
    list(APPEND AL4SAN_libs_to_find "al4san_starpu")
  elseif (QUARK_FOUND)
    list(APPEND AL4SAN_libs_to_find "al4san_quark")
      elseif (PARSEC_FOUND)
    list(APPEND AL4SAN_libs_to_find "al4san_parsec")
      elseif (OPENMP_FOUND)
    list(APPEND AL4SAN_libs_to_find "al4san_openmp")
  endif()
  list(APPEND AL4SAN_libs_to_find "coreblas")

  # call cmake macro to find the lib path
  if(AL4SAN_LIBDIR)
    foreach(al4san_lib ${AL4SAN_libs_to_find})
      set(AL4SAN_${al4san_lib}_LIBRARY "AL4SAN_${al4san_lib}_LIBRARY-NOTFOUND")
      find_library(AL4SAN_${al4san_lib}_LIBRARY
	NAMES ${al4san_lib}
	HINTS ${AL4SAN_LIBDIR})
    endforeach()
  else()
    if(AL4SAN_DIR)
      foreach(al4san_lib ${AL4SAN_libs_to_find})
	set(AL4SAN_${al4san_lib}_LIBRARY "AL4SAN_${al4san_lib}_LIBRARY-NOTFOUND")
	find_library(AL4SAN_${al4san_lib}_LIBRARY
	  NAMES ${al4san_lib}
	  HINTS ${AL4SAN_DIR}
	  PATH_SUFFIXES lib lib32 lib64)
      endforeach()
    else()
      foreach(al4san_lib ${AL4SAN_libs_to_find})
	set(AL4SAN_${al4san_lib}_LIBRARY "AL4SAN_${al4san_lib}_LIBRARY-NOTFOUND")
	find_library(AL4SAN_${al4san_lib}_LIBRARY
	  NAMES ${al4san_lib}
	  HINTS ${_lib_env})
      endforeach()
    endif()
  endif()

  # If found, add path to cmake variable
  # ------------------------------------
  foreach(al4san_lib ${AL4SAN_libs_to_find})

    get_filename_component(${al4san_lib}_lib_path ${AL4SAN_${al4san_lib}_LIBRARY} PATH)
    # set cmake variables (respects naming convention)
    if (AL4SAN_LIBRARIES)
      list(APPEND AL4SAN_LIBRARIES "${AL4SAN_${al4san_lib}_LIBRARY}")
    else()
      set(AL4SAN_LIBRARIES "${AL4SAN_${al4san_lib}_LIBRARY}")
    endif()
    if (AL4SAN_LIBRARY_DIRS)
      list(APPEND AL4SAN_LIBRARY_DIRS "${${al4san_lib}_lib_path}")
    else()
      set(AL4SAN_LIBRARY_DIRS "${${al4san_lib}_lib_path}")
    endif()
    mark_as_advanced(AL4SAN_${al4san_lib}_LIBRARY)

  endforeach(al4san_lib ${AL4SAN_libs_to_find})

  # check a function to validate the find
  if(AL4SAN_LIBRARIES)

    set(REQUIRED_LDFLAGS)
    set(REQUIRED_INCDIRS)
    set(REQUIRED_LIBDIRS)
    set(REQUIRED_LIBS)

    # AL4SAN
    if (AL4SAN_INCLUDE_DIRS)
      set(REQUIRED_INCDIRS "${AL4SAN_INCLUDE_DIRS}")
    endif()
    foreach(libdir ${AL4SAN_LIBRARY_DIRS})
      if (libdir)
	list(APPEND REQUIRED_LIBDIRS "${libdir}")
      endif()
    endforeach()
    set(REQUIRED_LIBS "${AL4SAN_LIBRARIES}")
    # STARPU
    if (STARPU_FOUND AND AL4SAN_LOOK_FOR_STARPU)
      if (STARPU_INCLUDE_DIRS_DEP)
	list(APPEND REQUIRED_INCDIRS "${STARPU_INCLUDE_DIRS_DEP}")
      elseif (STARPU_INCLUDE_DIRS)
	list(APPEND REQUIRED_INCDIRS "${STARPU_INCLUDE_DIRS}")
      endif()
      if(STARPU_LIBRARY_DIRS_DEP)
	list(APPEND REQUIRED_LIBDIRS "${STARPU_LIBRARY_DIRS_DEP}")
      elseif(STARPU_LIBRARY_DIRS)
	list(APPEND REQUIRED_LIBDIRS "${STARPU_LIBRARY_DIRS}")
      endif()
      if (STARPU_LIBRARIES_DEP)
	list(APPEND REQUIRED_LIBS "${STARPU_LIBRARIES_DEP}")
      elseif (STARPU_LIBRARIES)
	foreach(lib ${STARPU_LIBRARIES})
	  if (EXISTS ${lib} OR ${lib} MATCHES "^-")
	    list(APPEND REQUIRED_LIBS "${lib}")
	  else()
	    list(APPEND REQUIRED_LIBS "-l${lib}")
	  endif()
	endforeach()
      endif()
    endif()
    # QUARK
    if (QUARK_FOUND AND AL4SAN_LOOK_FOR_QUARK)
      if (QUARK_INCLUDE_DIRS_DEP)
	list(APPEND REQUIRED_INCDIRS "${QUARK_INCLUDE_DIRS_DEP}")
      elseif(QUARK_INCLUDE_DIRS)
	list(APPEND REQUIRED_INCDIRS "${QUARK_INCLUDE_DIRS}")
      endif()
      if(QUARK_LIBRARY_DIRS_DEP)
	list(APPEND REQUIRED_LIBDIRS "${QUARK_LIBRARY_DIRS_DEP}")
      elseif(QUARK_LIBRARY_DIRS)
	list(APPEND REQUIRED_LIBDIRS "${QUARK_LIBRARY_DIRS}")
      endif()
      if (QUARK_LIBRARY_DIRS_DEP)
	list(APPEND REQUIRED_LIBS "${QUARK_LIBRARIES_DEP}")
      elseif (QUARK_LIBRARY_DIRS_DEP)
	list(APPEND REQUIRED_LIBS "${QUARK_LIBRARIES}")
      endif()
    endif()
    # PARSEC
    if (PARSEC_FOUND AND AL4SAN_LOOK_FOR_PARSEC)
      if (PARSEC_INCLUDE_DIRS_DEP)
        list(APPEND REQUIRED_INCDIRS "${PARSEC_INCLUDE_DIRS_DEP}")
      elseif (PARSEC_INCLUDE_DIRS)
        list(APPEND REQUIRED_INCDIRS "${PARSEC_INCLUDE_DIRS}")
      endif()
      if(PARSEC_LIBRARY_DIRS_DEP)
        list(APPEND REQUIRED_LIBDIRS "${PARSEC_LIBRARY_DIRS_DEP}")
      elseif(PARSEC_LIBRARY_DIRS)
        list(APPEND REQUIRED_LIBDIRS "${PARSEC_LIBRARY_DIRS}")
      endif()
      if (PARSEC_LIBRARIES_DEP)
        list(APPEND REQUIRED_LIBS "${PARSEC_LIBRARIES_DEP}")
      elseif (PARSEC_LIBRARIES)
        foreach(lib ${PARSEC_LIBRARIES})
          if (EXISTS ${lib} OR ${lib} MATCHES "^-")
            list(APPEND REQUIRED_LIBS "${lib}")
          else()
            list(APPEND REQUIRED_LIBS "-l${lib}")
          endif()
        endforeach()
      endif()
    endif()
    # OPENMP
    if (OPENMP_FOUND AND AL4SAN_LOOK_FOR_OPENMP)
      if (OPENMP_INCLUDE_DIRS_DEP)
  list(APPEND REQUIRED_INCDIRS "${OPENMP_INCLUDE_DIRS_DEP}")
      elseif(OPENMP_INCLUDE_DIRS)
  list(APPEND REQUIRED_INCDIRS "${OPENMP_INCLUDE_DIRS}")
      endif()
      if(OPENMP_LIBRARY_DIRS_DEP)
  list(APPEND REQUIRED_LIBDIRS "${OPENMP_LIBRARY_DIRS_DEP}")
      elseif(OPENMP_LIBRARY_DIRS)
  list(APPEND REQUIRED_LIBDIRS "${OPENMP_LIBRARY_DIRS}")
      endif()
      if (OPENMP_LIBRARY_DIRS_DEP)
  list(APPEND REQUIRED_LIBS "${OPENMP_LIBRARIES_DEP}")
      elseif (OPENMP_LIBRARY_DIRS_DEP)
  list(APPEND REQUIRED_LIBS "${OPENMP_LIBRARIES}")
      endif()
    endif()    
    # CUDA
    if (CUDA_FOUND AND AL4SAN_LOOK_FOR_CUDA)
      if (CUDA_INCLUDE_DIRS)
	list(APPEND REQUIRED_INCDIRS "${CUDA_INCLUDE_DIRS}")
      endif()
      foreach(libdir ${CUDA_LIBRARY_DIRS})
	if (libdir)
	  list(APPEND REQUIRED_LIBDIRS "${libdir}")
	endif()
      endforeach()
      list(APPEND REQUIRED_LIBS "${CUDA_CUBLAS_LIBRARIES};${CUDA_LIBRARIES}")
    endif()
    # MAGMA
    if (MAGMA_FOUND AND AL4SAN_LOOK_FOR_MAGMA)
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
    endif()
    # MPI
    if (MPI_FOUND AND AL4SAN_LOOK_FOR_MPI)
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
    # HWLOC
    if (HWLOC_FOUND)
      if (HWLOC_INCLUDE_DIRS)
	list(APPEND REQUIRED_INCDIRS "${HWLOC_INCLUDE_DIRS}")
      endif()
      foreach(libdir ${HWLOC_LIBRARY_DIRS})
	if (libdir)
	  list(APPEND REQUIRED_LIBDIRS "${libdir}")
	endif()
      endforeach()
      foreach(lib ${HWLOC_LIBRARIES})
	if (EXISTS ${lib} OR ${lib} MATCHES "^-")
	  list(APPEND REQUIRED_LIBS "${lib}")
	else()
	  list(APPEND REQUIRED_LIBS "-l${lib}")
	endif()
      endforeach()
    endif()
    # TMG
    if (TMG_FOUND)
      if (TMG_INCLUDE_DIRS_DEP)
	list(APPEND REQUIRED_INCDIRS "${TMG_INCLUDE_DIRS_DEP}")
      elseif (TMG_INCLUDE_DIRS)
	list(APPEND REQUIRED_INCDIRS "${TMG_INCLUDE_DIRS}")
      endif()
      if(TMG_LIBRARY_DIRS_DEP)
	list(APPEND REQUIRED_LIBDIRS "${TMG_LIBRARY_DIRS_DEP}")
      elseif(TMG_LIBRARY_DIRS)
	list(APPEND REQUIRED_LIBDIRS "${TMG_LIBRARY_DIRS}")
      endif()
      if (TMG_LIBRARIES_DEP)
	list(APPEND REQUIRED_LIBS "${TMG_LIBRARIES_DEP}")
      elseif(TMG_LIBRARIES)
	list(APPEND REQUIRED_LIBS "${TMG_LIBRARIES}")
      endif()
      if (TMG_LINKER_FLAGS)
	list(APPEND REQUIRED_LDFLAGS "${TMG_LINKER_FLAGS}")
      endif()
    endif()
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
    list(APPEND REQUIRED_LIBS ${AL4SAN_EXTRA_LIBRARIES})

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
    unset(AL4SAN_WORKS CACHE)
    include(CheckFunctionExists)
    check_function_exists(AL4SAN_Init AL4SAN_WORKS)
    mark_as_advanced(AL4SAN_WORKS)

    if(AL4SAN_WORKS)
      # save link with dependencies
      set(AL4SAN_LIBRARIES_DEP "${REQUIRED_LIBS}")
      set(AL4SAN_LIBRARY_DIRS_DEP "${REQUIRED_LIBDIRS}")
      set(AL4SAN_INCLUDE_DIRS_DEP "${REQUIRED_INCDIRS}")
      set(AL4SAN_LINKER_FLAGS "${REQUIRED_LDFLAGS}")
      list(REMOVE_DUPLICATES AL4SAN_LIBRARY_DIRS_DEP)
      list(REMOVE_DUPLICATES AL4SAN_INCLUDE_DIRS_DEP)
      list(REMOVE_DUPLICATES AL4SAN_LINKER_FLAGS)
    else()
      if(NOT AL4SAN_FIND_QUIETLY)
	message(STATUS "Looking for al4san : test of AL4SAN_Init fails")
	message(STATUS "CMAKE_REQUIRED_LIBRARIES: ${CMAKE_REQUIRED_LIBRARIES}")
	message(STATUS "CMAKE_REQUIRED_INCLUDES: ${CMAKE_REQUIRED_INCLUDES}")
	message(STATUS "Check in CMakeFiles/CMakeError.log to figure out why it fails")
	message(STATUS "Maybe AL4SAN is linked with specific libraries. "
	  "Have you tried with COMPONENTS (STARPU/QUARK, CUDA, MAGMA, MPI, FXT)? "
	  "See the explanation in FindAL4SAN.cmake.")
      endif()
    endif()
    set(CMAKE_REQUIRED_INCLUDES)
    set(CMAKE_REQUIRED_FLAGS)
    set(CMAKE_REQUIRED_LIBRARIES)
  endif(AL4SAN_LIBRARIES)

endif( (NOT PKG_CONFIG_EXECUTABLE) OR (PKG_CONFIG_EXECUTABLE AND NOT AL4SAN_FOUND) OR (AL4SAN_GIVEN_BY_USER) )

if (AL4SAN_LIBRARIES)
  if (AL4SAN_LIBRARY_DIRS)
    foreach(dir ${AL4SAN_LIBRARY_DIRS})
      if ("${dir}" MATCHES "al4san")
	set(first_lib_path "${dir}")
      endif()
    endforeach()
  else()
    list(GET AL4SAN_LIBRARIES 0 first_lib)
    get_filename_component(first_lib_path "${first_lib}" PATH)
  endif()
  if (${first_lib_path} MATCHES "/lib(32|64)?$")
    string(REGEX REPLACE "/lib(32|64)?$" "" not_cached_dir "${first_lib_path}")
    set(AL4SAN_DIR_FOUND "${not_cached_dir}" CACHE PATH "Installation directory of AL4SAN library" FORCE)
  else()
    set(AL4SAN_DIR_FOUND "${first_lib_path}" CACHE PATH "Installation directory of AL4SAN library" FORCE)
  endif()
endif()
mark_as_advanced(AL4SAN_DIR)
mark_as_advanced(AL4SAN_DIR_FOUND)

# check that AL4SAN has been found
# ---------------------------------
include(FindPackageHandleStandardArgs)
if (PKG_CONFIG_EXECUTABLE AND AL4SAN_FOUND)
  find_package_handle_standard_args(AL4SAN DEFAULT_MSG
    AL4SAN_LIBRARIES)
else()
  find_package_handle_standard_args(AL4SAN DEFAULT_MSG
    AL4SAN_LIBRARIES
    AL4SAN_WORKS)
endif()
