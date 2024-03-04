
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file ImportHCore.cmake
# @brief Checks for the Hcore library and includes it in the project if it is not already present.
# @version 1.1.0
# @author Mahmoud ElKarargy
# @author Sameh Abdulah
# @date 2023-03-15

#Configurations
set(name "HCORE")
set(tag "v0.1.3")
set(version "0.1.3")
set(flag "")
set(is_cmake ON)
set(is_git ON)
set(auto_gen OFF)
set(url "https://github.com/ecrc/hcore.git")

include(macros/ImportDependency)
ImportDependency(${name} ${tag} ${version} ${url} "${flag}" "" ${is_cmake} ${is_git} ${auto_gen})

message(STATUS "${name} done")
