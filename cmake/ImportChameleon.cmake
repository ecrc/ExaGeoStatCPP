
# Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file ImportChameleon.cmake
# @brief This script checks for Chameleon and includes it in the project if it is not already a target.
# @version 1.1.0
# @author Mahmoud ElKarargy
# @author Sameh Abdulah
# @date 2023-03-13

# Set basic configuration variables for the Chameleon library.
# 'name' is set to "CHAMELEON", which is the identifier used for the dependency throughout the script.
set(name "CHAMELEON")
# 'tag' specifies the version tag of the Chameleon library to fetch, indicating a specific state of the source code in the repository.
set(tag "v1.1.0")
# 'version' specifies the version of the Chameleon library. This may be used to ensure compatibility or meet specific requirements.
set(version "1.1.0")
# 'flag' can be used to pass additional flags to the configure/make commands when building the dependency. Flags here are for CUDA, BLAS vendor, MPI, StarPU, and disabling testing.
set(flag -DCHAMELEON_USE_CUDA=${USE_CUDA} \-DBLA_VENDOR=${BLA_VENDOR} \-DCHAMELEON_USE_MPI=${USE_MPI} \-DCHAMELEON_SCHED_STARPU=ON \-DCHAMELEON_ENABLE_TESTING=OFF)
# 'is_cmake' is a boolean flag indicating whether the Chameleon library uses CMake for its build system. It's set to ON, meaning it does.
set(is_cmake ON)
# 'is_git' is a boolean flag indicating whether the Chameleon library's source code is hosted in a git repository. It's set to ON.
set(is_git ON)
# 'auto_gen' is a boolean flag indicating whether to use autogen scripts for the configuration process. It's set to OFF here.
set(auto_gen OFF)
# 'url' specifies the location of the Chameleon library's source code repository.
set(url "https://gitlab.inria.fr/solverstack/chameleon.git")

# Include the 'ImportDependency' macro script located in the 'macros' directory.
# This macro is responsible for importing and possibly installing the dependency.
include(macros/ImportDependency)
# Call the 'ImportDependency' macro with the previously set configuration parameters.
# This macro checks if CHAMELEON is already available; if not, it proceeds to fetch, configure, build, and install it.
ImportDependency(${name} ${tag} ${version} ${url} "${flag}" "" ${is_cmake} ${is_git} ${auto_gen})

# Additional include directories for Chameleon are specified, ensuring the project can find Chameleon's headers.
# The AFTER keyword specifies that these directories should be searched after the default ones.
include_directories(AFTER ${CHAMELEON_DIR_FOUND}/include/coreblas)

if(NOT CHAMELEON_DIR_FOUND)
    set(CHAMELEON_DIR_FOUND ${CMAKE_INSTALL_PREFIX}/CHAMELEON)
else()
    include_directories(${CHAMELEON_DIR_FOUND})
endif()
message(STATUS "CHAMELEON_DIR_FOUND: ${CHAMELEON_DIR_FOUND}")
include_directories(${CHAMELEON_DIR_FOUND}/chameleon-src)
# Print a status message indicating the completion of Chameleon's inclusion process.
message(STATUS "${name} done")