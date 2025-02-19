
# Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file FindHICMA-X.cmake
# @brief This is a CMakeLists file for finding HiCMA-X and link and include it's headers
# @version 2.0.0
# @author Mahmoud ElKarargy
# @date 2024-09-28

# Include pkg-config
find_package(PkgConfig QUIET)

# Try to find dplasma and parsec via pkg-config
if(PKG_CONFIG_FOUND)
    pkg_check_modules(DPLASMA_PKG dplasma)
    pkg_check_modules(PARSEC_PKG parsec)
    if(DPLASMA_PKG_FOUND AND PARSEC_PKG_FOUND)
        # Try to find the HICMA-X or hicma-x path in the library directories
        string(FIND "${PARSEC_PKG_LIBRARY_DIRS}" "HICMA-X" HICMA_X_START)
        if(HICMA_X_START EQUAL -1)
            string(FIND "${PARSEC_PKG_LIBRARY_DIRS}" "hicma-x" HICMA_X_START)
        endif()
        if(HICMA_X_START GREATER -1)
            # Extract the full path to HICMA-X or hicma-x and set the include directory
            string(REGEX MATCH "([^;]*(HICMA-X|hicma-x)[^;]*/lib)" HICMA_X_LIB_PATH "${PARSEC_PKG_LIBRARY_DIRS}")
            get_filename_component(HICMA_X_ROOT "${HICMA_X_LIB_PATH}" DIRECTORY) # Go one level up
            set(HICMA-X_INCLUDE_DIRS "${HICMA_X_ROOT}/include")  # Set the include path
        endif()
        # TODO: This is not generalized for the case of hicma installed manually
        set(HICMA_X_SRC_DIR ${HICMA_X_ROOT}/hicma-x-src)
        set(HICMA_X_BIN_DIR ${HICMA_X_ROOT}/bin)
        set(HICMA-X_FOUND TRUE)
        set(HICMA-X_LIBRARIES ${DPLASMA_PKG_LIBRARIES} ${PARSEC_PKG_LIBRARIES})
        set(HICMA-X_LIBRARY_DIRS "${HICMA_X_LIB_PATH}")
        # Add a search for lib64 directories and set HICMA-X_LIBRARY_DIRS_DEP
        set(HICMA-X_LIBRARY_DIRS_DEP "${HICMA_X_LIB_PATH}64")

        find_library(HICMA_PARSEC_LIB hicma_parsec PATHS ${HICMA-X_LIBRARY_DIRS_DEP})

        if(HICMA_PARSEC_LIB)
            list(APPEND HICMA-X_LIBRARIES ${HICMA_PARSEC_LIB})
        else()
            message(FATAL_ERROR "libhicma_parsec.so not found")
        endif()

    endif()
endif()

# Fallback: Manual search if pkg-config fails or HICMA-X path isn't set
if(NOT HICMA-X_FOUND)
    # Improved search to handle multiple possible paths and fallback for include directories
    find_path(HICMA-X_INCLUDE_DIR
            NAMES hicma.h
            PATHS
            ${CMAKE_CURRENT_LIST_DIR}/../hicma-x/include
            /usr/local/include/hicma-x
            /usr/local/include
            /usr/include/hicma-x
            /usr/include
            DOC "Path to HICMA-X include directory"
            )

    # Search for the main HICMA-X library
    find_library(HICMA-X_LIBRARY
            NAMES hicma-x
            PATHS
            ${CMAKE_CURRENT_LIST_DIR}/../hicma-x/lib
            /usr/local/lib
            /usr/lib
            DOC "Path to HICMA-X library"
            )

    # Search for the hicma_parsec library in the lib64 directory if it's not found in the standard lib
    find_library(HICMA_PARSEC_LIB
            NAMES hicma_parsec
            PATHS
            ${CMAKE_CURRENT_LIST_DIR}/../hicma-x/lib64
            /usr/local/lib64
            /usr/lib64
            DOC "Path to HICMA-Parsec library"
            )

    # Check if both the include directory and libraries were found
    if(HICMA-X_INCLUDE_DIR AND HICMA-X_LIBRARY AND HICMA_PARSEC_LIB)
        set(HICMA-X_FOUND TRUE)
        # Combine the found libraries
        set(HICMA-X_LIBRARIES ${HICMA-X_LIBRARY} ${HICMA_PARSEC_LIB})
        # Set the include directory
        set(HICMA-X_INCLUDE_DIRS "${HICMA-X_INCLUDE_DIR}")
        # Include both lib and lib64 directories
        # TODO: This paths are not generalized, if the install is not with the same dir.
        set(HICMA-X_LIBRARY_DIRS "${HICMA-X_LIBRARY}/lib")
        set(HICMA-X_LIBRARY_DIRS_DEP "${HICMA-X_LIBRARY}/lib64")
    else()
        set(HICMA-X_FOUND FALSE)
    endif()
endif()

# Mark the variables as advanced to keep the CMake GUI clean
mark_as_advanced(HICMA-X_INCLUDE_DIR HICMA-X_LIBRARY HICMA_PARSEC_LIB)

# Provide feedback on whether the library was found
if(HICMA-X_FOUND)
    message(STATUS "Found HICMA-X")
else()
    message("Could not find HICMA-X or its dependencies (dplasma, parsec)")
endif()
