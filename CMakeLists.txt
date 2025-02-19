
# Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @brief This is a CMakeLists.txt file for a C++ project ExaGeoStat.
# The project is a parallel high performance unified framework for geographical statistics on manycore systems.
# The file sets up variables and finds dependencies required for the project.
# It also provides options to enable building tests, building examples, building documentation, and enabling a packaging system for distribution.
# @version 1.1.0
# @author Mahmoud ElKarargy
# @author Sameh Abdulah
# @date 2024-02-04

# Set the minimum CMake version required to 3.20
cmake_minimum_required(VERSION 3.20 FATAL_ERROR)
cmake_policy(SET CMP0048 NEW)

# Set project options
option(USE_CUDA "Use CUDA, if available" OFF)
option(USE_MPI "Use MPI, if available" OFF)
option(BUILD_TESTS "Option to enable building tests" OFF)
option(BUILD_HEAVY_TESTS "Option to enable building heavy tests, This may take a lot of time" OFF)
option(BUILD_EXAMPLES "Option to enable building examples" ON)
option(BUILD_DOCS "Build documentation in docs directory" ON)
option(USE_R "Enable the use of R and Rcpp in the project" OFF)
option(CREATE_PACKAGE "Enable a packaging system for distribution" OFF)

# Declare the RUNTIME_TYPE variable with a default value
set(RUNTIME_TYPE "starpu" CACHE STRING "Specify the runtime type (e.g., starpu or parsec)")
# Cmake Module Paths
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake;${CMAKE_MODULE_PATH}")
if (${BUILD_SHARED_LIBS})
    set(BLA_STATIC OFF)
else ()
    set(BLA_STATIC ON)
endif ()

# Project Name and Version
project(ExaGeoStatCPP VERSION 2.0.0 DESCRIPTION "ExaGeoStatCPP is a parallel high performance unified framework for geostatistics on manycore systems.")
set(CMAKE_CXX_EXTENSIONS OFF)
string(TOUPPER ${RUNTIME_TYPE} RUNTIME_TYPE)

if(RUNTIME_TYPE STREQUAL "PARSEC")
    message(STATUS "MPI is required to be enabled in order to use PaRSEC")
    set(USE_MPI ON)
    set(USE_HICMA OFF)
endif()

# Select toolchain based on whether CUDA is enabled or not
if (USE_CUDA)
    message("")
    message("---------------------------------------- CUDA")
    # Enable CUDA and include CudaToolchain
    add_definitions(-DUSE_CUDA=TRUE)
    enable_language(CUDA)
    include(toolchains/CudaToolchain)
    # Set BLA_VENDOR to NVHPC for CUDA-enabled builds
    set(BLA_VENDOR NVHPC)
    list(APPEND STARPU_COMPONENT_LIST "CUDA")
else ()
    message("-- Build x86 Support")
    # Include GccToolchain for non-CUDA builds - Gcc
    include(toolchains/GccToolchain)
endif ()

# Show the current version of CMake.
message(STATUS "CMAKE VERSION: ${CMAKE_VERSION}")
# Enable C++ language
enable_language(CXX)
# Get the current path of the project.
add_compile_definitions(PROJECT_SOURCE_DIR="${PROJECT_SOURCE_DIR}/")
# Set kernels path.
add_definitions(-DKERNELS_PATH="${PROJECT_SOURCE_DIR}/inst/include/kernels/concrete/")

# ExaGeoStatCPP depends on CUDA
# -------------------------------
if (USE_CUDA)
    message("-- Build CUDA Support")
else ()
    message("-- Build x86 Support")
    set(gpu_backend CACHE STRING "none" FORCE)
    unset(BLA_VENDOR)
endif ()

# ExaGeoStatCPP depends on MPI
# -------------------------------
if (USE_MPI)
    message("")
    message("---------------------------------------- MPI")
    # Enable MPI and include MPI
    add_definitions(-DUSE_MPI=TRUE)
    message(STATUS "Trying to find MPI")
    find_package(MPI REQUIRED)
    include_directories(${MPI_INCLUDE_PATH})
    list(APPEND LIBS ${MPI_LIBRARIES})
    list(APPEND STARPU_COMPONENT_LIST "MPI")
endif ()

# ExaGeoStatCPP depends on LAPACKE
#-----------------------------
message("")
message("---------------------------------------- LAPACKE")
find_package(LAPACKE)
list(APPEND LIBS ${LAPACKE_LIBRARIES})
link_directories(${LAPACKE_LIBRARY_DIRS_DEP})
include_directories(${LAPACKE_INCLUDE_DIRS})

if (${BLA_VENDOR} MATCHES "Intel10_64lp" OR ${BLA_VENDOR} MATCHES "Intel" OR ${BLA_VENDOR} MATCHES "Intel10_64lp_seq")
    message(STATUS "Linked to MKL")
    add_compile_definitions(USE_MKL)
endif()

# Add all dependencies for ExaGeoStatCPP
#-----------------------------

# Print installation path of ExaGeoStatCPP.
message(STATUS "Installation path : ${CMAKE_INSTALL_PREFIX}")

# ExaGeoStatCPP depends on HWLoc
# -------------------------------
include(ImportHwloc)
list(APPEND STARPU_COMPONENT_LIST "HWLOC")
string(REPLACE ";" " " STARPU_COMPONENT_STRING "${STARPU_COMPONENT_LIST}")

# ExaGeoStatCPP depends on NLOPT
# -------------------------------
include(ImportNLOPT)

# ExaGeoStatCPP depends on GSL
# -------------------------------
include(ImportGSL)

message("---------------------------------------- ${RUNTIME_TYPE}")
if(RUNTIME_TYPE STREQUAL "STARPU")
    message(STATUS "Using StarPU as the runtime")
    # ExaGeoStatCPP depends on StarPU runtime
    # -------------------------------
    include(ImportStarPu)
    # ExaGeoStatCPP depends on HiCMA
    # -------------------------------
    if (USE_HICMA)
        add_definitions(-DUSE_HICMA=TRUE)
        include(ImportHCore)
        include(ImportStarsH)
        include(ImportHiCMA)
    endif ()
    # ExaGeoStatCPP depends on Chameleon
    # -------------------------------
    include(ImportChameleon)

elseif(RUNTIME_TYPE STREQUAL "PARSEC")
    message(STATUS "Using PaRSEC as the runtime")
    # ExaGeoStatCPP depends on StarPU runtime
    # -------------------------------
    include(ImportStarsH)
    include(ImportHCore)
    include(ImportHiCMAX)
else()
    message(FATAL_ERROR "Unknown RUNTIME_TYPE: ${RUNTIME_TYPE}. Supported values are 'STARPU' or 'PARSEC'.")
endif()

# ExaGeoStatCPP depends on LAPACK/BLASPP
# -------------------------------
include(ImportBLASPP)
include(ImportLapack)

# ExaGeoStatCPP Documentation
if (BUILD_DOCS)
    find_package(Doxygen)
    if (DOXYGEN_FOUND)
        add_subdirectory("docs")
    else ()
        message(STATUS "Doxygen NOT found, skipping it")
    endif ()
endif ()

# Include directories for ExaGeoStatCPP
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/inst/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/prerequisites)
set(MY_LOGGER_PATH ${CMAKE_CURRENT_SOURCE_DIR})
add_definitions(-DMY_LOGGER_PATH="${CMAKE_CURRENT_SOURCE_DIR}")

if (USE_R)
    message("")
    message("---------------------------------------- Rcpp")
    # Find R and Rcpp using FindR Module
    find_package(R REQUIRED)
    if (${R_FOUND})
        message(STATUS "Using R technology")
        list(APPEND LIBS R)
        add_definitions(-DUSE_R)
    endif ()
endif ()

# Add src Directory to expose added libraries
add_subdirectory(src)

# Define USE_STARPU or USE_RUNTIME based on the selected runtime
if(RUNTIME_TYPE STREQUAL "STARPU")
    target_compile_definitions(${PROJECT_NAME} PUBLIC DEFAULT_RUNTIME=1)
elseif(RUNTIME_TYPE STREQUAL "PARSEC")
    target_compile_definitions(${PROJECT_NAME} PUBLIC DEFAULT_RUNTIME=0)
endif()

# Creates a new INTERFACE library target named ${PROJECT_NAME}_INTERFACE.
# The INTERFACE keyword specifies that this library will not be built, but instead will only be used for its properties.
add_library(${PROJECT_NAME}_INTERFACE INTERFACE)
target_link_libraries(${PROJECT_NAME}_INTERFACE INTERFACE ${PROJECT_NAME})

# Add linker options to the target
target_link_options(${PROJECT_NAME}_INTERFACE INTERFACE "SHELL:-Wl,--whole-archive $<TARGET_FILE:${PROJECT_NAME}> -Wl,--no-whole-archive")

# Add tests if enabled
if (${BUILD_TESTS})
    message(STATUS "Building Tests")
    include(ImportCatch2)
    include(Catch)
    include(CTest)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp-tests)
    enable_testing()
endif ()

# Add heavy tests if enabled
if (${BUILD_HEAVY_TESTS})
    message(STATUS "Building Heavy Tests")
    include(ImportCatch2)
    include(Catch)
    include(CTest)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/tests/heavy-tests)
    enable_testing()
    message(STATUS "Building heavy tests will enable examples too")
    set(BUILD_EXAMPLES ON)
endif ()

# Add examples if enabled
if (BUILD_EXAMPLES)
    message(STATUS "Building Examples is Enabled")
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/examples)
endif ()

# ExaGeoStatCPP Dependence export messages for users
# -------------------------------
message(" \n \t ** Configurations of ExaGeoStatCPP and installation of dependence is done successfully ** ")
message("\t - Export the following line to avoid re-install dependencies each time. -")
message("\t ----------------------------------------------------------------------------------------------------------------------------------- ")
message("\t export PKG_CONFIG_PATH=${CMAKE_INSTALL_PREFIX}/HWLOC/lib/pkgconfig:$PKG_CONFIG_PATH")
message("\t export PKG_CONFIG_PATH=${CMAKE_INSTALL_PREFIX}/NLOPT/lib/pkgconfig:${CMAKE_INSTALL_PREFIX}/NLOPT/lib64/pkgconfig:$PKG_CONFIG_PATH")
message("\t export PKG_CONFIG_PATH=${CMAKE_INSTALL_PREFIX}/GSL/lib/pkgconfig:$PKG_CONFIG_PATH")
if(RUNTIME_TYPE STREQUAL "STARPU")
    message("\t export PKG_CONFIG_PATH=${CMAKE_INSTALL_PREFIX}/STARPU/lib/pkgconfig:$PKG_CONFIG_PATH")
    message("\t export PKG_CONFIG_PATH=${CMAKE_INSTALL_PREFIX}/CHAMELEON/lib/pkgconfig:$PKG_CONFIG_PATH")
    if(USE_HICMA)
        message("\t export PKG_CONFIG_PATH=${CMAKE_INSTALL_PREFIX}/STARSH/lib/pkgconfig:$PKG_CONFIG_PATH")
        message("\t export PKG_CONFIG_PATH=${CMAKE_INSTALL_PREFIX}/HCORE/lib/pkgconfig:$PKG_CONFIG_PATH")
        message("\t export PKG_CONFIG_PATH=${CMAKE_INSTALL_PREFIX}/HICMA/lib/pkgconfig:$PKG_CONFIG_PATH")
    endif()
elseif(RUNTIME_TYPE STREQUAL "PARSEC")
    message("\t export PKG_CONFIG_PATH=${CMAKE_INSTALL_PREFIX}/STARSH/lib/pkgconfig:$PKG_CONFIG_PATH")
    message("\t export PKG_CONFIG_PATH=${CMAKE_INSTALL_PREFIX}/HCORE/lib/pkgconfig:$PKG_CONFIG_PATH")
    message("\t export PKG_CONFIG_PATH=${CMAKE_INSTALL_PREFIX}/HICMA-X/lib/pkgconfig:${CMAKE_INSTALL_PREFIX}/HICMA-X/lib64:$PKG_CONFIG_PATH")
endif()
message("\t ----------------------------------------------------------------------------------------------------------------------------------- \n")

# Installation of ExaGeoStatCPP
install(DIRECTORY inst/include/ DESTINATION ${CMAKE_INSTALL_PREFIX}/EXAGEOSTATCPP/include)

## Install cmake find package.
include(CMakePackageConfigHelpers)
write_basic_package_version_file("${CMAKE_INSTALL_PREFIX}/EXAGEOSTATCPP/${PROJECT_NAME}ConfigVersion.cmake" COMPATIBILITY ExactVersion)
install(
        FILES
        "${CMAKE_INSTALL_PREFIX}/EXAGEOSTATCPP/${PROJECT_NAME}ConfigVersion.cmake"
        DESTINATION ${CMAKE_INSTALL_PREFIX}/EXAGEOSTATCPP/lib/cmake/${PROJECT_NAME}
)

configure_file(${PROJECT_NAME}Config.cmake.in
        ${CMAKE_INSTALL_PREFIX}/EXAGEOSTATCPP/${PROJECT_NAME}Config.cmake @ONLY)

install(
        FILES
        "${CMAKE_INSTALL_PREFIX}/EXAGEOSTATCPP/${PROJECT_NAME}Config.cmake"
        DESTINATION ${CMAKE_INSTALL_PREFIX}/EXAGEOSTATCPP/lib/cmake/${PROJECT_NAME}
)

install(
        DIRECTORY
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake"
        DESTINATION ${CMAKE_INSTALL_PREFIX}/EXAGEOSTATCPP/Modules
)

## Generate pkg-config file
configure_file(package.pc.in
        ${CMAKE_INSTALL_PREFIX}/EXAGEOSTATCPP/lib/pkgconfig/${PROJECT_NAME}.pc @ONLY)
install(
        FILES
        "${CMAKE_INSTALL_PREFIX}/EXAGEOSTATCPP/lib/pkgconfig/${PROJECT_NAME}.pc"
        DESTINATION ${CMAKE_INSTALL_PREFIX}/EXAGEOSTATCPP/lib/pkgconfig/
)

if (CREATE_PACKAGE)
    ##################
    # Release source #
    ##################
    set(CPACK_SOURCE_GENERATOR "TGZ")
    set(CPACK_PACKAGE_NAME "${PROJECT_NAME}")
    set(CPACK_PACKAGE_DESCRIPTION_FILE ${CMAKE_CURRENT_SOURCE_DIR}/README.md)
    set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "ExaGeoStatCPP is a parallel high performance unified framework for geostatistics on manycore systems. Its abbreviation stands for 'Exascale Geostatistics'.")
    set(CPACK_PACKAGE_VERSION "${${PROJECT_NAME}_VERSION}")
    set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}")
    set(CPACK_SOURCE_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}")
    set(CPACK_PACKAGE_VENDOR "KAUST")
    set(CPACK_PACKAGE_CONTACT "sameh.abdulah@kaust.edu.sa")
    set(CPACK_RESOURCE_FILE_README ${CMAKE_CURRENT_SOURCE_DIR}/README.md)
    set(CPACK_RESOURCE_FILE_LICENSE ${CMAKE_CURRENT_SOURCE_DIR}/LICENSE)
    set(CPACK_SOURCE_IGNORE_FILES "bin;.git;Jenkinsfile")
    include(CPack)

    message("\t - Export the following line if you want to use ExaGeoStatCPP in another software. -")
    message("\t ----------------------------------------------------------------------------------------------------------------------------------------------------- ")
    message("\t export PKG_CONFIG_PATH=${CMAKE_INSTALL_PREFIX}/EXAGEOSTATCPP/lib/pkgconfig:$PKG_CONFIG_PATH")
    message("\t ----------------------------------------------------------------------------------------------------------------------------------------------------- ")
endif ()