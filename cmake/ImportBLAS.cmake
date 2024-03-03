# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file ImportBLAS.cmake
# @brief This file searches for the BLAS library and includes it if not already included.
# @version 1.1.0
# @author Mahmoud ElKarargy
# @date 2023-03-12

#Configurations
set(name "BLAS")
set(tag "v0.3.21")
set(version "0.3.21")
set(flag "")
set(is_cmake ON)
set(is_git ON)
set(auto_gen OFF)
set(url "https://github.com/xianyi/OpenBLAS")

include(macros/ImportDependency)
ImportDependency(${name} ${tag} ${version} ${url} "${flag}" "" ${is_cmake} ${is_git} ${auto_gen})
message(STATUS "${name} done")
