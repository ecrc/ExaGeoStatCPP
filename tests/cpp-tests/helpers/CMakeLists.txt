# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @version 1.0.0
# @author Mahmoud ElKarargy
# @author Sameh Abdulah
# @date 2023-12-08

set(EXAGEOSTAT_TESTFILES

        ${CMAKE_CURRENT_SOURCE_DIR}/TestPredictionHelpers.cpp

        ${EXAGEOSTAT_TESTFILES}
        PARENT_SCOPE
        )
