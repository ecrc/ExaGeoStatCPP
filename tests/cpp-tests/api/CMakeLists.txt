
# Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @version 1.1.0
# @author Mahmoud ElKarargy
# @author Sameh Abdulah
# @date 2023-08-08

set(EXAGEOSTAT_TESTFILES

        ${CMAKE_CURRENT_SOURCE_DIR}/TestExaGeoStatApi.cpp

        ${EXAGEOSTAT_TESTFILES}
        PARENT_SCOPE
        )
