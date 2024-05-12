# Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @version 1.1.0
# @author Mahmoud ElKarargy
# @date 2024-02-09

set(EXAGEOSTAT_TESTFILES

        ${CMAKE_CURRENT_SOURCE_DIR}/TestAllRFunctions.cpp
        ${EXAGEOSTAT_TESTFILES}
        PARENT_SCOPE
        )

