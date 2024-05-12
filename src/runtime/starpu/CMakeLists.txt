
# Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @version 1.1.0
# @brief CMake build script for StarPu runtime
# @author Mahmoud ElKarargy
# @date 2024-02-19

# Include the concrete implementations of the StarPu codelets classes
file(GLOB ALL_CODELETS ${CMAKE_CURRENT_SOURCE_DIR}/concrete/*.cpp)

# Include StarPu helpers.
add_subdirectory(helpers)

# Define the sources for the library
set(SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/StarPuFunctions.cpp
        ${ALL_CODELETS}
        ${SOURCES}
        PARENT_SCOPE
        )

# Automatically add new codelets header files to StarPuCodeletsHeaders.hpp

# File documentation
set(DOCUMENTATION_STRING "
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file StarPuCodelets.hpp
 * @brief Header file to include all codelet classes.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @date 2024-02-25
**/

")

set(PATH_TO_CODELETS_HEADERS "${PROJECT_SOURCE_DIR}/inst/include/runtime/starpu")
set(GENERATED_INCLUDES "${PATH_TO_CODELETS_HEADERS}/StarPuCodeletsHeaders.hpp")

file(GLOB_RECURSE HEADER_FILES "${PATH_TO_CODELETS_HEADERS}/concrete/*.hpp")
file(WRITE ${GENERATED_INCLUDES} ${DOCUMENTATION_STRING})

foreach (HEADER_FILE ${HEADER_FILES})
    # Construct the include directive by stripping the known base directory part from the full path
    string(REPLACE "${PATH_TO_CODELETS_HEADERS}" "" HEADER_RELATIVE_PATH ${HEADER_FILE})
    # Use angle brackets and the desired base path for the include directive
    file(APPEND ${GENERATED_INCLUDES} "#include <runtime/starpu${HEADER_RELATIVE_PATH}>\n")
endforeach ()
