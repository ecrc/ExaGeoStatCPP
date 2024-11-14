
# Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @brief Find and include STARSH library as a dependency.
# @version 2.0.0
# @author Mahmoud ElKarargy
# @author Sameh Abdulah
# @date 2024-09-28

# Configuration parameters for integrating the STARSH library
# 'name' is set to "STARSH" to identify the STARSH library within this script.
set(name "STARSH")

# Check the value of RUNTIME_TYPE and configure STARSH accordingly
if(RUNTIME_TYPE STREQUAL "STARPU")
    # Default values for STARPU runtime
    set(STARSH_TAG "v0.3.1")
    set(STARSH_VERSION "0.3.1")
    set(STARSH_URL "https://github.com/ecrc/stars-h.git")
    message(STATUS "RUNTIME_TYPE is STARPU. Using default STARSH configuration.")

elseif(RUNTIME_TYPE STREQUAL "PARSEC")
    # Custom values for PARSEC runtime
    set(STARSH_TAG "sabdulah/non-gaussian-kernel")
    set(STARSH_VERSION "0")
    set(STARSH_URL "https://github.com/SAbdulah/stars-h.git")
    message(STATUS "RUNTIME_TYPE is PARSEC. Using custom STARSH configuration for PARSEC.")
endif()

# 'flag' is used for additional build configuration options, specifically disabling StarPU and optionally enabling MPI.
set(flag \-DSTARPU=OFF \-DMPI=${USE_MPI} \-DBLA_VENDOR=${BLA_VENDOR})
# 'is_cmake' indicates that STARSH uses CMake as its build system, set to ON.
set(is_cmake ON)
# 'is_git' denotes that the source code for STARSH is hosted on a Git repository, set to ON.
set(is_git ON)
# 'auto_gen' signals whether autogen scripts are needed for the build process; it is set to OFF for STARSH.
set(auto_gen OFF)
# The 'ImportDependency' macro, located in the 'macros' directory, is included to manage the import and setup of the STARSH library.
include(macros/ImportDependency)
# The 'ImportDependency' macro is called with the configuration parameters set above to manage the detection, fetching, and setup of STARSH.
ImportDependency(${name} ${STARSH_TAG} ${STARSH_VERSION} ${STARSH_URL} "${flag}" "" ${is_cmake} ${is_git} ${auto_gen})

# A message is output to indicate the successful integration of the STARSH library into the project.
message(STATUS "${name} done")
