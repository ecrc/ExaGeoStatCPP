
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file ImportChameleon.cmake
# @brief This script checks for Chameleon and includes it in the project if it is not already a target.
# @version 1.0.1
# @author Mahmoud ElKarargy
# @date 2023-03-13

#Configurations
set(name "Catch2")
set(tag "v3.3.2")
set(version "3.3.2")
set(flag "")
set(is_cmake ON)
set(is_git ON)
set(auto_gen OFF)
set(url "https://github.com/catchorg/Catch2.git")

include(macros/ImportDependency)
ImportDependency(${name} ${tag} ${version} ${url} "${flag}" "" ${is_cmake} ${is_git} ${auto_gen})
message(STATUS "${name} done")
