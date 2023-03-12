
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# Copyright (c) 2023 by Brightskies inc,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @version 1.0.0
# @author Sameh Abdulah
# @date 2023-01-31

set(EXAGEOSTAT_TESTFILES ${CMAKE_CURRENT_SOURCE_DIR}/test_main.cpp)

add_subdirectory(configurations/data-generation)
add_subdirectory(data-generators)

enable_testing()

add_executable(exageostat-tests ${EXAGEOSTAT_TESTFILES})
target_link_libraries(exageostat-tests ${PROJECT_NAME})
