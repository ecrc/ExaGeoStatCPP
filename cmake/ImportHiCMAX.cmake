# Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file ImportHiCMAX.cmake
# @brief Find and include HiCMA-X library as a dependency.
# @version 2.0.0
# @author Mahmoud ElKarargy
# @author Sameh Abdulah
# @date 2024-09-21

# Configuration settings for integrating the HICMA-X library into the project
# 'name' sets the identifier for the HICMA-X library within this script to "HICMA-X".
set(name "HICMA-X")
# Set the version tag for HiCMA-X.
set(tag "FIX-package-installation-MK")
# Flags to configure the build for HiCMA-X, including precision settings for DPLASMA
# and disabling GPU support for both CUDA and HIP.
set(flags '-DDPLASMA_PRECISIONS="s;d"' \-DPARSEC_WITH_DEVEL_HEADERS=ON \-DCMAKE_Fortran_FLAGS="-Wno-main"
    \-DPARSEC_GPU_WITH_HIP=OFF \-DPARSEC_GPU_WITH_CUDA=OFF \-DPARSEC_HAVE_CUDA=OFF \-DPARSEC_DIST_SHORT_LIMIT=0
    \-DPARSEC_DIST_COLLECTIVES=ON \-DPARSEC_HAVE_DEV_CUDA_SUPPORT=OFF \-DDPLASMA_HAVE_CUDA=OFF \-DBLA_VENDOR=${BLA_VENDOR})
# Indicates that HiCMA-X uses CMake for its build system.
set(is_cmake ON)
# Indicates that HiCMA-X is hosted on a Git repository.
set(is_git ON)
# Indicates that autogen scripts are not required for HiCMA-X.
set(auto_gen OFF)
# Set the URL of the HiCMA-X GitHub repository.
set(url "https://github.com/SAbdulah/hicma-x-dev.git")
# Include the macro to import HiCMA-X as a dependency.
include(macros/ImportDependency)

# Use the ImportDependency macro to handle fetching, detecting, and setting up HiCMA-X.
ImportDependency(${name} ${tag} "" ${url} "${flags}" "" ${is_cmake} ${is_git} ${auto_gen})

# Include necessary directories for HiCMA-X and its dependencies.
include_directories(${HICMA_X_SRC_DIR})
include_directories(${HICMA_X_SRC_DIR}/dplasma/src)
include_directories(${HICMA_X_SRC_DIR}/hicma_parsec)
include_directories(${HICMA_X_SRC_DIR}/bin/dplasma/src)
# Display a status message indicating that HiCMA-X has been successfully included.
message(STATUS "HiCMA-X done")
