
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @version 1.0.0
# @author Sameh Abdulah
# @date 2023-01-30

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/configurations)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/data-generators)

set(LIB_NAME ${PROJECT_NAME})
add_library(${LIB_NAME}
        STATIC
        ${SOURCES}
        )
target_link_libraries(${LIB_NAME} ${LIBS})
target_compile_definitions(${LIB_NAME} PUBLIC ${COMPILE_DEFINITIONS})

## Installation steps for headers and libraries.
set_target_properties(${LIB_NAME}
        PROPERTIES
        VERSION ${${PROJECT_NAME}_VERSION}
        )
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
