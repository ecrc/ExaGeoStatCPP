
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @version 1.0.0
# @author Mahmoud ElKarargy
# @date 2023-01-31

add_subdirectory(api)
add_subdirectory(linear-algebra-solvers)
add_subdirectory(data-generators)
add_subdirectory(configurations/data-generation)
add_subdirectory(kernels)
add_subdirectory(helpers)

enable_testing()
add_executable(exageostat-tests ${EXAGEOSTAT_TESTFILES})
target_link_libraries(exageostat-tests Catch2::Catch2WithMain ${PROJECT_NAME}_INTERFACE)
catch_discover_tests(exageostat-tests)