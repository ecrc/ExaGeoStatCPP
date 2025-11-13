
# Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file ImportPnetCDF.cmake
# @brief Checks for the PnetCDF library and includes it in the project if it is not already present.
# @version 2.0.0
# @author Mahmoud ElKarargy
# @author Sameh Abdulah
# @date 2024-11-14

# Compute MPI installation root from the MPI C++ compiler path
if (DEFINED MPI_CXX_COMPILER)
  get_filename_component(MPI_COMPILER_DIR "${MPI_CXX_COMPILER}" DIRECTORY)        # .../bin
  get_filename_component(MPI_ROOT "${MPI_COMPILER_DIR}" DIRECTORY)               # parent of bin
else()
  message(FATAL_ERROR "MPI_CXX_COMPILER not set; cannot locate MPI installation")
endif()

# Configuration settings for the integration of the NLOPT library
# 'name' is assigned to "NLOPT", serving as the identifier for this library within the script.
set(name "PnetCDF")
# 'tag' defines "PnetCDF-1_12_0" as the version tag of NLOPT, indicating the specific release to be utilized.
set(tag "tag.v1.10.0")
# 'version' specifies "1.12.0" as the version of the NLOPT library, ensuring compatibility with the project's requirements.
set(version "1.10.0")
# 'flag' is intended for additional configuration options during the build process. A space is placed as a placeholder.
set(flag \--enable-shared \--with-mpi=${MPI_ROOT})
# 'is_cmake' indicates that NLOPT uses CMake for its build system, which is set to ON.
set(is_cmake OFF)
# 'is_git' denotes that the NLOPT source code is hosted in a Git repository, which is set to ON.
set(is_git ON)
# 'auto_gen' signals whether autogen scripts are required for the build process, which is set to OFF for NLOPT.
set(auto_gen ON)
# 'url' provides the location of the NLOPT source code repository on GitHub.
set(url "https://github.com/Parallel-NetCDF/PnetCDF")

# The 'ImportDependency' macro script, located in the 'macros' directory, is included for managing the import and setup of the NLOPT library.
include(macros/ImportDependency)
# The 'ImportDependency' macro is invoked with the above-defined parameters to handle the detection, fetching, and integration of NLOPT into the project.
ImportDependency(${name} ${tag} ${version} ${url} "${flag}" "" ${is_cmake} ${is_git} ${auto_gen})

# A status message is outputted to indicate the successful integration of the NLOPT library into the project.
message(STATUS "${name} done")
