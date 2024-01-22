
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @brief Find and include STARSH library as a dependency.
# @version 1.0.1
# @author Mahmoud ElKarargy
# @author Sameh Abdulah
# @date 2023-03-13

#Configurations
set(name "STARSH")
set(tag "v0.3.1")
set(version "0.3.1")
set(flag \-DSTARPU=OFF \-DMPI=${USE_MPI})
set(is_cmake ON)
set(is_git ON)
set(auto_gen OFF)
set(url "https://github.com/ecrc/stars-h.git")

include(macros/ImportDependency)
ImportDependency(${name} ${tag} ${version} ${url} "${flag}" "" ${is_cmake} ${is_git} ${auto_gen})

message(STATUS "${name} done")
