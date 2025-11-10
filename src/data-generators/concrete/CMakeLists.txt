
# Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @brief CMake configuration file for the Synthetic Generator module
# @version 1.1.0
# @author Mahmoud ElKarargy
# @author Sameh Abdulah
# @date 2023-02-14

# Add source files to the parent scope - always include SyntheticGenerator
list(APPEND SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/SyntheticGenerator.cpp)

if (USE_CLIMATE_EMULATOR)
    if (RUNTIME_TYPE STREQUAL "STARPU")
        # Use StarPU/CHAMELEON version
        list(APPEND SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/StageZeroGenerator.cpp)
    elseif (RUNTIME_TYPE STREQUAL "PARSEC")
        # Use PaRSEC version
        list(APPEND SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/StageZeroGeneratorParsec.cpp)
        list(APPEND SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/runtime/ParsecGenerator.cpp)
    endif()
endif()

# Set the final SOURCES to parent scope
set(SOURCES ${SOURCES} PARENT_SCOPE)
