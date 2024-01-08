# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file ImportBlasPP.cmake
# @brief This file searches for the BLAS++ library and includes it if not already included.
# @version 1.0.1
# @author Mahmoud ElKarargy
# @date 2023-03-12

include(ImportBlas)

#Configurations
set(name blaspp)
string(TOUPPER ${name} capital_name)
set(tag "v2023.01.00")
# Set installation flags
if (USE_CUDA)
    set(flag "-Dgpu_backend=cuda")
else()
    set(flag "-Dgpu_backend=")
endif ()

set(version "2023.01.00")
set(is_cmake ON)
set(is_git ON)
set(auto_gen OFF)
set(url "https://github.com/icl-utk-edu/blaspp")

set(${name}_DIR "${CMAKE_INSTALL_PREFIX}/${capital_name}/lib/cmake/${name}")
include(macros/ImportDependency)
ImportDependency(${name} ${tag} ${version} ${url} "${flag}" ${is_cmake} ${is_git} ${auto_gen})
message(STATUS "${name} done")
