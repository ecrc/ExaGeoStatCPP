
# Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @version 1.1.0
# @author Mahmoud ElKarargy
# @date 2023-12-20

enable_testing()
add_executable(exageostat-heavy-tests ${CMAKE_CURRENT_SOURCE_DIR}/HeavyTests.cpp ${CMAKE_CURRENT_SOURCE_DIR}/ExamplesTests.cpp)
target_link_libraries(exageostat-heavy-tests Catch2::Catch2WithMain ${PROJECT_NAME}_INTERFACE)
catch_discover_tests(exageostat-heavy-tests)