
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file ImportHwloc.cmake
# @brief Find and include Hwloc library as a dependency.
# @version 1.0.1
# @author Mahmoud ElKarargy
# @author Sameh Abdulah
# @date 2023-03-15

#Configurations
set(name "HWLOC")
set(tag "hwloc-2.4.0")
set(version "2.4.0")
set(flag "")
set(is_cmake OFF)
set(is_git ON)
set(auto_gen ON)
set(url "https://github.com/open-mpi/hwloc")

include(macros/ImportDependency)
ImportDependency(${name} ${tag} ${version} ${url} "${flag}" ${is_cmake} ${is_git} ${auto_gen})

message(STATUS "${name} done")
