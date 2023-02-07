# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @version 1.0.0
# @author Sameh Abdulah
# @date 2023-01-31

# REQUIRED FOR TESTS TO LINK LIBRARIES
if (NOT EXISTS "${PROJECT_SOURCE_DIR}/libs/ecrc_cmake/modules")
    find_package(Git REQUIRED)
    execute_process(COMMAND ${GIT_EXECUTABLE} submodule init WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} RESULT_VARIABLE _res_init OUTPUT_QUIET ERROR_QUIET)
    execute_process(COMMAND ${GIT_EXECUTABLE} submodule update WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} RESULT_VARIABLE _res_update OUTPUT_QUIET ERROR_QUIET)
    if (${_res_init} GREATER 0 OR ${_res_update} GREATER 0)
        message(FATAL_ERROR "ECRC CMake modules were not found.\n"
                "We tried: 'git submodule init && git submodule update' and resulted in error")
    endif ()
endif ()

# ECRC INITIALIZATION
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/libs/ecrc_cmake/modules)
set(ECRC_CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules/ecrc/modules)

include(EcrcInit)