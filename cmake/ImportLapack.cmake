
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file ImportLapack.cmake
# @brief Find and include LAPACK library as a dependency.
# @version 1.0.1
# @author Mahmoud ElKarargy
# @author Sameh Abdulah
# @date 2023-03-12

#Configurations
set(name "LAPACK")
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
