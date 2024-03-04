
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @version 1.1.0
# @author Mahmoud ElKarargy
# @date 2024-01-24

add_subdirectory(api)
add_subdirectory(configurations)
add_subdirectory(data-generators)
add_subdirectory(hardware)
add_subdirectory(helpers)
add_subdirectory(kernels)
add_subdirectory(linear-algebra-solvers)
add_subdirectory(prediction)
add_subdirectory(results)

if(USE_HICMA)
    add_subdirectory(data-units)
endif ()

if (USE_R)
    add_subdirectory(Rcpp-adapters)
endif ()

enable_testing()
add_executable(exageostat-tests ${EXAGEOSTAT_TESTFILES})
target_link_libraries(exageostat-tests Catch2::Catch2WithMain ${PROJECT_NAME}_INTERFACE)
catch_discover_tests(exageostat-tests)