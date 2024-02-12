
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @version 1.1.0
# @author Mahmoud ElKarargy
# @date 2023-04-29

# Automatically add all kernels tests in the concrete directory.
file(GLOB ALL_KERNELS ${CMAKE_CURRENT_SOURCE_DIR}/concrete/*.cpp)

set(EXAGEOSTAT_TESTFILES
        ${ALL_KERNELS}
        ${EXAGEOSTAT_TESTFILES}
        PARENT_SCOPE
        )