
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file ImportChameleon.cmake
# @brief This script checks for Chameleon and includes it in the project if it is not already a target.
# @version 1.1.0
# @author Mahmoud ElKarargy
# @author Sameh Abdulah
# @date 2023-03-13

#Configurations
set(name "CHAMELEON")
set(tag "v1.1.0")
set(version "1.1.0")
set(flag -DCHAMELEON_USE_CUDA=${USE_CUDA} \-DBLA_VENDOR=${BLA_VENDOR} \-DCHAMELEON_USE_MPI=${USE_MPI} \-DCHAMELEON_SCHED_STARPU=ON \-DCHAMELEON_ENABLE_TESTING=OFF)
set(is_cmake ON)
set(is_git ON)
set(auto_gen OFF)
set(url "https://gitlab.inria.fr/solverstack/chameleon.git")

include(macros/ImportDependency)
ImportDependency(${name} ${tag} ${version} ${url} "${flag}" "" ${is_cmake} ${is_git} ${auto_gen})

include_directories(AFTER ${CHAMELEON_DIR_FOUND}/include/coreblas)
include_directories(${CHAMELEON_DIR_FOUND}/chameleon-src)

message(STATUS "${name} done")
