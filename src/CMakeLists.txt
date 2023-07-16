
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @brief This file contains the CMake configuration for the ExaGeoStat library.
# @version 1.0.0
# @author Sameh Abdulah
# @date 2023-01-30

# Add subdirectories for configurations, data-generators, data-units, and linear-algebra-solvers.
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/configurations)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/data-generators)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/data-units)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/linear-algebra-solvers)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/kernels)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/api)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/helpers)

# Set the name of the library to be created.
set(LIB_NAME ${PROJECT_NAME})

# Create the library with the specified source files and linking libraries.
add_library(${LIB_NAME}
        STATIC
        ${SOURCES}
        )
target_compile_definitions(${LIB_NAME} PUBLIC ${COMPILE_DEFINITIONS})
target_link_libraries(${LIB_NAME} PUBLIC ${LIBS})

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
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        )
install(EXPORT ${LIB_NAME}CoreConfig DESTINATION lib/cmake/${PROJECT_NAME})