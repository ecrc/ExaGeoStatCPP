
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @brief This is a CMakeLists.txt file for a C++ project ExaGeoStat.
# The project is a parallel high performance unified framework for geographical statistics on many core systems.
# The file sets up variables and finds dependencies required for the project.
# It also provides options to enable building tests, building examples, building documentation, and enabling a packaging system for distribution.
# @version 1.0.0
# @author Sameh Abdulah
# @date 2023-01-30

# Set the minimum CMake version required to 3.20
cmake_minimum_required(VERSION 3.20 FATAL_ERROR)
cmake_policy(SET CMP0048 NEW)

# Cmake Module Paths
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake;${CMAKE_MODULE_PATH}")

# Select toolchain based on whether CUDA is enabled or not
if (USE_CUDA)
    # Enable CUDA and include CudaToolchain
    add_definitions(-DUSE_CUDA=TRUE)
    enable_language(CUDA)
    include(toolchains/CudaToolchain)
    # Set BLA_VENDOR to NVHPC for CUDA-enabled builds
    set(BLA_VENDOR NVHPC)
else ()
    message("-- Build x86 Support")
    # Include GccToolchain for non-CUDA builds - Gcc
    include(toolchains/GccToolchain)
endif ()

# Project Name and Version
project(exageostatcpp VERSION 1.0.0 DESCRIPTION "ExaGeoStat is a parallel high performance unified framework for geostatistics on manycore systems.")

# Set project options
option(USE_CUDA "Use Cuda, if available" false)
option(USE_MPI "Use MPI, if available" false)
option(EXAGEOSTAT_BUILD_TESTS "Option to enable building tests" ON)
option(EXAGEOSTAT_BUILD_EXAMPLES "Option to enable building examples" ON)
option(EXAGEOSTAT_BUILD_DOCS "Build documentation in docs directory" ON)
option(EXAGEOSTAT_PACKAGE "Enable a packaging system for distribution" OFF)

# Show the current version of CMake.
message(STATUS "CMAKE VERSION: ${CMAKE_VERSION}")
# Enable C++ language
enable_language(CXX)

# Add all dependencies for ExaGeoStat PP
include(ImportGFortran)
include(ImportOpenMP)
if (USE_CUDA)
    message("-- Build CUDA Support")
    # Import CuSolver and Lapack libraries
    include(ImportCuSolver)
    include(ImportBlasPP)
    unset(BLA_VENDOR)
    include(ImportLapack)

    set(LIBS
            LAPACK::LAPACK
            ${LIBS}
            )
else ()
    message("-- Build x86 Support")
    set(gpu_backend CACHE "none" FORCE)
    unset(BLA_VENDOR)
endif ()

# Add all dependencies for Exageostat
# EXAGEOSTAT depends on a MPI
# -------------------------------
if (USE_MPI)
    message(STATUS "Trying to find MPI")
    find_package(MPI REQUIRED)
    list(APPEND STARPU_COMPONENT_LIST "MPI")
endif ()

# EXAGEOSTAT depends on a Hwloc
# -------------------------------
include(ImportHwloc)
set(STARPU_COMPONENT_LIST "HWLOC")

# EXAGEOSTAT depends on a runtime
# -------------------------------
include(ImportStarPu)

# EXAGEOSTAT depends on a GSL
# -------------------------------
include(ImportGSL)

# EXAGEOSTAT depends on a NLOPT
# -------------------------------
include(ImportNLOPT)

# EXAGEOSTAT depends on HiCMA
# -------------------------------
if (EXAGEOSTAT_USE_HiCMA)
    add_definitions(-DEXAGEOSTAT_USE_HiCMA=TRUE)
    message(STATUS "Add Hcore, Dependency needed for HiCMA")
    include(ImportHcore)
    message(STATUS "Add StarsH, Dependency needed for HiCMA")
    include(ImportStarsH)
    include(ImportHiCMA)
endif ()

# EXAGEOSTAT depends on CHAMELEON
# -------------------------------
if (EXAGEOSTAT_USE_CHAMELEON)
    add_definitions(-DEXAGEOSTAT_USE_CHAMELEON=TRUE)
    include(ImportChameleon)
endif ()

# EXAGEOSTAT depends on a LAPACKE
# -------------------------------
include(ImportBlas)
include(ImportLapack)

find_package(LAPACKE)
list(APPEND LIBS ${LAPACKE_LIBRARIES})
link_directories(${LAPACKE_LIBRARY_DIRS_DEP})
include_directories(${LAPACKE_INCLUDE_DIRS})


# EXAGEOSTAT DOCUMENTATIONS
if (EXAGEOSTAT_BUILD_DOCS)
    find_package(Doxygen)

    if (DOXYGEN_FOUND)
        add_subdirectory("docs")
    else ()
        message(STATUS "Doxygen NOT found, skipping it")
    endif ()
endif ()

# Include directories for Exageostat-cpp
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/inst/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/prerequisites)
set(MY_LOGGER_PATH ${CMAKE_CURRENT_SOURCE_DIR})
add_definitions(-DMY_LOGGER_PATH="${CMAKE_CURRENT_SOURCE_DIR}")

# Add src Directory to expose added libraries
add_subdirectory(src)

# Creates a new INTERFACE library target named ${PROJECT_NAME}_INTERFACE.
# The INTERFACE keyword specifies that this library will not be built, but instead will only be used for its properties.
add_library(${PROJECT_NAME}_INTERFACE INTERFACE)
target_link_libraries(${PROJECT_NAME}_INTERFACE INTERFACE ${PROJECT_NAME})
# Add linker options to the target
target_link_options(${PROJECT_NAME}_INTERFACE INTERFACE "SHELL:-Wl,--whole-archive $<TARGET_FILE:${PROJECT_NAME}> -Wl,--no-whole-archive")

# Install headers
install(TARGETS exageostatcpp
        DESTINATION lib/
        PUBLIC_HEADER DESTINATION include/
        )
# Add tests if enabled
if (${EXAGEOSTAT_BUILD_TESTS})
    message(STATUS "Building Tests is Enabled")
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp-tests)
    include(CTest)
    add_test(NAME exageostat-tests
            COMMAND exageostat-tests
            )
endif ()

if (EXAGEOSTAT_BUILD_EXAMPLES)
    message(STATUS "Building Examples  is Enabled")
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/examples)
endif ()

# Print installation path of Exageostat.
message(STATUS "Installation path : ${EXAGEOSTAT_INSTALL_PREFIX}")
# Installation actions
install(DIRECTORY include/${PROJECT_NAME} DESTINATION include)
## Install cmake find package.
include(CMakePackageConfigHelpers)
write_basic_package_version_file("${CMAKE_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake" COMPATIBILITY ExactVersion)
install(
        FILES
        "${CMAKE_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
        DESTINATION lib/cmake/${PROJECT_NAME}
)

configure_file(${PROJECT_NAME}Config.cmake.in
        ${PROJECT_NAME}Config.cmake @ONLY)

install(
        FILES
        "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
        DESTINATION lib/cmake/${PROJECT_NAME}
)

install(
        DIRECTORY
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake"
        DESTINATION lib/cmake/${PROJECT_NAME}/Modules
)

## Generate pkg-config file
configure_file(package.pc.in
        lib/pkgconfig/${PROJECT_NAME}.pc @ONLY)
install(
        FILES
        "${PROJECT_BINARY_DIR}/lib/pkgconfig/${PROJECT_NAME}.pc"
        DESTINATION lib/pkgconfig/
)

if (EXAGEOSTAT_PACKAGE)
    ##################
    # Release source #
    ##################
    set(CPACK_SOURCE_GENERATOR "TGZ")
    set(CPACK_PACKAGE_NAME "${PROJECT_NAME}")
    set(CPACK_PACKAGE_DESCRIPTION_FILE ${CMAKE_CURRENT_SOURCE_DIR}/README.md)
    set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "ExaGeoStat is a parallel high performance unified framework for geostatistics on manycore systems. Its abbreviation stands for 'Exascale Geostatistics'.")
    set(CPACK_PACKAGE_VERSION "${${PROJECT_NAME}_VERSION}")
    set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}")
    set(CPACK_SOURCE_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}")
    set(CPACK_PACKAGE_VENDOR "KAUST")
    set(CPACK_PACKAGE_CONTACT "sameh.abdulah@kaust.edu.sa")
    set(CPACK_RESOURCE_FILE_README ${CMAKE_CURRENT_SOURCE_DIR}/README.md)
    set(CPACK_RESOURCE_FILE_LICENSE ${CMAKE_CURRENT_SOURCE_DIR}/LICENSE)
    set(CPACK_SOURCE_IGNORE_FILES "bin;.git;.gitmodules;Jenkinsfile")
    include(CPack)
endif ()

message(" \n \t ** Configurations of ExaGeoStat and installation of dependence is done successfully ** ")
message("\t Export the following line to avoid re-install dependencies each time, This line assume that you haven't changed the installation path. ")
message("\t If not, Please change the following paths with your installation path. \n")
message("\t ------------------------------------------------------------------------------------------------------------------------------- ")
message("\t export PKG_CONFIG_PATH=$PWD/installdir/_deps/CHAMELEON/lib/pkgconfig:$PWD/installdir/_deps/GSL/lib/pkgconfig:$PKG_CONFIG_PATH")
message("\t export PKG_CONFIG_PATH=$PWD/installdir/_deps/HCORE/lib/pkgconfig:$PWD/installdir/_deps/HICMA/lib/pkgconfig:$PKG_CONFIG_PATH")
message("\t export PKG_CONFIG_PATH=$PWD/installdir/_deps/HWLOC/lib/pkgconfig:$PWD/installdir/_deps/NLOPT/lib/pkgconfig:$PKG_CONFIG_PATH")
message("\t export PKG_CONFIG_PATH=$PWD/installdir/_deps/STARPU/lib/pkgconfig:$PWD/installdir/_deps/STARSH/lib/pkgconfig:$PKG_CONFIG_PATH")
message("\t ------------------------------------------------------------------------------------------------------------------------------- \n")
