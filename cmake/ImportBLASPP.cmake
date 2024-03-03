# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file ImportBLASPP.cmake
# @brief This file searches for the BLAS++ library and includes it if not already included.
# @version 1.0.1
# @author Mahmoud ElKarargy
# @date 2023-03-12

include(ImportBLAS)

#Configurations
set(name blaspp)
set(tag "v2023.01.00")
set(url "https://github.com/icl-utk-edu/blaspp")
set(${name}_DIR "${CMAKE_INSTALL_PREFIX}/${capital_name}/lib/cmake/${name}")

include(FetchContent)
FetchContent_Declare(${name} GIT_REPOSITORY "${url}" GIT_TAG "${tag}")
FetchContent_MakeAvailable(${name})

set(LIBS ${name} ${LIBS})
message(STATUS "${name} done")