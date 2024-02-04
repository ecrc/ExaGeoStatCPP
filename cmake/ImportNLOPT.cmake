
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file ImportNLOPT.cmake
# @brief Find and include NLOPT library as a dependency.
# @version 1.1.0
# @author Mahmoud ElKarargy
# @author Sameh Abdulah
# @date 2023-03-26

#Configurations
set(name "NLOPT")
set(tag "v2.7.1")
set(version "2.7.1")
set(flag " ")
set(is_cmake ON)
set(is_git ON)
set(auto_gen OFF)
set(url "https://github.com/stevengj/nlopt")

include(macros/ImportDependency)
ImportDependency(${name} ${tag} ${version} ${url} "${flag}" "" ${is_cmake} ${is_git} ${auto_gen})

message(STATUS "${name} done")
