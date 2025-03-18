
# Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @version 1.1.0
# @brief CMake build script for Parsec directory.
# @author Mahmoud ElKarargy
# @date 2024-03-10

# Add subdirectory for JDF files
add_subdirectory(jdf)

# Set source files for ClimateEmulator
set(SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/ParsecFunctions.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/JDFHelperFunctions.c
        ${SOURCES}
        PARENT_SCOPE
)

# Set JDF generated source files for Climate Emulator
set(JDF_GENERATED_SOURCES
        ${JDF_GENERATED_SOURCES}
        PARENT_SCOPE
)

