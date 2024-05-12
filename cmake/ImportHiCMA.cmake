
# Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file ImportHiCMA.cmake
# @brief Find and include HiCMA library as a dependency.
# @version 1.1.0
# @author Mahmoud ElKarargy
# @author Sameh Abdulah
# @date 2023-03-13

# Configuration settings for integrating the HICMA library into the project
# 'name' sets the identifier for the HICMA library within this script to "HICMA".
set(name "HICMA")
# 'tag' specifies "v1.0.0" as the version tag of HICMA, denoting a specific release to be fetched.
set(tag "v1.0.0")
# 'version' defines "1.0.0" as the version of the HICMA library, ensuring compatibility with project requirements.
set(version "1.0.0")
# 'flag' is used to pass additional configuration options during the build, specifically for enabling or disabling MPI support based on the project's needs.
set(flag -DHICMA_USE_MPI=${USE_MPI})
# 'is_cmake' indicates that HICMA uses CMake for its build system, set to ON.
set(is_cmake ON)
# 'is_git' denotes that HICMA's source code is hosted on a Git repository, set to ON.
set(is_git ON)
# 'auto_gen' signals whether autogen scripts are necessary for the build process; it is set to OFF for HICMA.
set(auto_gen OFF)
# 'url' provides the GitHub repository URL for HICMA, specifying where the source code can be cloned from.
set(url "https://github.com/ecrc/hicma.git")

# The 'ImportDependency' macro, located in the 'macros' directory, is included. This macro handles the import and setup of the HICMA library.
include(macros/ImportDependency)
# The 'ImportDependency' macro is invoked with the previously defined parameters to manage the detection, fetching, and setup of HICMA.
ImportDependency(${name} ${tag} ${version} ${url} "${flag}" "" ${is_cmake} ${is_git} ${auto_gen})

# Directories containing HiCMA headers are included in the project, ensuring that HiCMA's functions and types are accessible.
include_directories(${HICMA_LIBDIR}/../hicma-src/hicma_ext)
include_directories(${HICMA_LIBDIR}/../hicma_ext)

# A status message is displayed to indicate the successful inclusion of the HiCMA library into the project.
message(STATUS "HiCMA done")

