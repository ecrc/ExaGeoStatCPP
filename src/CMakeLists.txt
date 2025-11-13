
# Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @brief This file contains the CMake configuration for the ExaGeoStat library.
# @version 1.1.0
# @author Mahmoud ElKarargy
# @author Sameh Abdulah
# @date 2024-02-04

# Add subdirectories for configurations, data-generators, data-units, and linear-algebra-solvers.
add_subdirectory(configurations)
add_subdirectory(kernels)
add_subdirectory(helpers)
add_subdirectory(data-units)
add_subdirectory(runtime)
add_subdirectory(api)
add_subdirectory(results)
add_subdirectory(hardware)
add_subdirectory(data-generators)
add_subdirectory(data-loader)
add_subdirectory(prediction)
add_subdirectory(runtime-solver)

# Check the value of RUNTIME_TYPE and configure chameleon/starpu accordingly
if (RUNTIME_TYPE STREQUAL "PARSEC")
    add_subdirectory(data-analyzer)
    add_subdirectory(data-transformer)
else()
    add_subdirectory(linear-algebra-solvers)
endif ()

if (USE_R)
    add_subdirectory(Rcpp-adapters)
endif ()

# Set the name of the library to be created.
set(LIB_NAME ${PROJECT_NAME})
# Create the library with the specified source files and linking libraries.
add_library(${LIB_NAME} SHARED ${SOURCES} ${JDF_GENERATED_SOURCES})

target_compile_definitions(${LIB_NAME} PUBLIC ${COMPILE_DEFINITIONS})
target_link_libraries(${LIB_NAME} ${LIBS})

# Set the version of the library.
set_target_properties(${LIB_NAME}
        PROPERTIES
        VERSION ${${PROJECT_NAME}_VERSION}
        )

# Set the installation directories for the library.
target_include_directories(${LIB_NAME}
        PUBLIC
        "$<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>"
        )
install(TARGETS ${LIB_NAME}
        EXPORT ${LIB_NAME}CoreConfig
        LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/EXAGEOSTATCPP/lib
        ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/EXAGEOSTATCPP/lib
        )

install(EXPORT ${LIB_NAME}CoreConfig DESTINATION ${CMAKE_INSTALL_PREFIX}/EXAGEOSTATCPP/lib/cmake/)