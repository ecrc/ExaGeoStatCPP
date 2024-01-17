
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file ImportHiCMA.cmake
# @brief Find and include HiCMA library as a dependency.
# @version 1.0.1
# @author Mahmoud ElKarargy
# @author Sameh Abdulah
# @date 2023-03-13

#Configurations
set(name "HICMA")
set(tag "v1.0.0")
set(version "1.0.0")
set(flag -DHICMA_USE_MPI=${USE_MPI} \-DCMAKE_C_FLAGS=-fPIC)
set(is_cmake ON)
set(is_git ON)
set(auto_gen OFF)
set(url "https://github.com/ecrc/hicma.git")

include(macros/ImportDependency)
ImportDependency(${name} ${tag} ${version} ${url} "${flag}" "" ${is_cmake} ${is_git} ${auto_gen})

# Include HiCMA headers in the project.
include_directories(${HICMA_LIBDIR}/../hicma-src/hicma_ext)
message(STATUS "HiCMA done")
