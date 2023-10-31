
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file CMakeLists.txt
# @brief This is a CMakeLists.txt file for a C++ project ExaGeoStat.
# The project is a parallel high performance unified framework for geographical statistics on manycore systems.
# The file sets up variables and finds dependencies required for the project.
# It also provides options to enable building tests, building examples, building documentation, and enabling a packaging system for distribution.
# @version 1.0.0
# @author Mahmoud ElKarargy
# @author Sameh Abdulah
# @date 2023-01-30

# Set the minimum CMake version required to 3.20
cmake_minimum_required(VERSION 3.20 FATAL_ERROR)
cmake_policy(SET CMP0048 NEW)

# Set project options
option(USE_CUDA "Use Cuda, if available" false)
option(USE_MPI "Use MPI, if available" false)
option(EXAGEOSTAT_BUILD_TESTS "Option to enable building tests" ON)
option(EXAGEOSTAT_BUILD_EXAMPLES "Option to enable building examples" ON)
option(EXAGEOSTAT_BUILD_DOCS "Build documentation in docs directory" ON)
option(EXAGEOSTAT_PACKAGE "Enable a packaging system for distribution" OFF)

# Cmake Module Paths
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake;${CMAKE_MODULE_PATH}")

# Project Name and Version
project(exageostatcpp VERSION 1.0.0 DESCRIPTION "ExaGeoStat is a parallel high performance unified framework for geostatistics on manycore systems.")

# Show the current version of CMake.
message(STATUS "CMAKE VERSION: ${CMAKE_VERSION}")

# Select toolchain based on whether CUDA is enabled or not
if (USE_CUDA)
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


# Enable C++ language
enable_language(CXX)

add_definitions(
	        -DLOG_PATH="${PROJECT_SOURCE_DIR}/synthetic_ds/"
		        -DKERNELS_PATH="${PROJECT_SOURCE_DIR}/inst/include/kernels/concrete/"
			)



# EXAGEOSTAT depends on CBLAS
#---------------------------
find_package(CBLAS COMPONENTS BLASEXT)
if(BLAS_FOUND)
	    if (BLAS_LIBRARY_DIRS)
		            # the RPATH to be used when installing
			            list(APPEND CMAKE_INSTALL_RPATH "${BLAS_LIBRARY_DIRS}")
				        endif()
					    if(BLAS_LINKER_FLAGS)
						            list(APPEND CMAKE_EXE_LINKER_FLAGS "${BLAS_LINKER_FLAGS}")
							        endif()
							else()
								    message(FATAL_ERROR "BLAS library has not been found")
							    endif()
							    if(CBLAS_FOUND)
								        include_directories(${CBLAS_INCLUDE_DIRS})
									    if(CBLAS_LIBRARY_DIRS)
										            # the RPATH to be used when installing
											            list(APPEND CMAKE_INSTALL_RPATH "${CBLAS_LIBRARY_DIRS}")
												        endif()
													    if (CBLAS_LIBRARIES)
														            if (CBLAS_LIBRARIES_DEP)
																                list(INSERT EXAGEOSTAT_DEP 0 ${CBLAS_LIBRARIES_DEP})
																		        else()
																				            list(INSERT EXAGEOSTAT_DEP 0 ${CBLAS_LIBRARIES})
																					            endif()
																						        endif()
																						else()
																							    if(ECRC_VERBOSE_FIND_PACKAGE)
																								            if(CBLAS_STANDALONE OR NOT CBLAS_WORKS)
																										                if (NOT CBLAS_cblas.h_DIRS)
																													                Print_Find_Header_Status(cblas cblas.h)
																															            endif ()
																																                if (NOT CBLAS_cblas_LIBRARY)
																																			                Print_Find_Library_Status(cblas libcblas)
																																					            endif ()
																																						            endif()
																																							        else()
																																									        message(WARNING "CBLAS library has not been found and ECRC_VERBOSE_FIND_PACKAGE is set to OFF."
																																											        " Try to activate ECRC_VERBOSE_FIND_PACKAGE option (-DECRC_VERBOSE_FIND_PACKAGE=ON) to get some hints for the detection")
																																											    endif()
																																											        message(FATAL_ERROR "A CBLAS library is required but has not been found")
																																											endif()


																																											find_package(LAPACKE)
																																											list(APPEND LIBS ${LAPACKE_LIBRARIES})
																																											link_directories(${LAPACKE_LIBRARY_DIRS_DEP})
																																											include_directories(${LAPACKE_INCLUDE_DIRS})


# Add all dependencies for ExaGeoStat PP
#include(ImportGFortran)
#include(ImportOpenMP)
if (USE_CUDA)
    message("-- Build CUDA Support")
    # Import CuSolver and Lapack libraries
    #    include(ImportCuSolver)
         include(ImportBlasPP)
     unset(BLA_VENDOR)
    include(ImportLapack)

    set(LIBS
           LAPACK::LAPACK
           ${LIBS}
           )
else ()
    message("-- Build x86 Support")
    set(gpu_backend CACHE STRING "none" FORCE)
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

# Check if no path is set for installation
if(NOT EXAGEOSTAT_INSTALL_PREFIX)
    message(FATAL_ERROR "Installation path not set! Please use -DEXAGEOSTAT_INSTALL_PREFIX=path/to/install or use ./config.sh")
endif()
# Print installation path of Exageostat.
message(STATUS "Installation path : ${EXAGEOSTAT_INSTALL_PREFIX}")

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
if (EXAGEOSTAT_USE_HICMA)
    add_definitions(-DEXAGEOSTAT_USE_HICMA=TRUE)
    message(STATUS "Add Hcore, Dependency needed for HiCMA")
    include(ImportHcore)
    message(STATUS "Add StarsH, Dependency needed for HiCMA")
    include(ImportStarsH)
    include(ImportHiCMA)
endif ()

# EXAGEOSTAT depends on CHAMELEON
# -------------------------------
include(ImportChameleon)

# EXAGEOSTAT depends on a LAPACKE
# -------------------------------
include(ImportBlas)
include(ImportLapack)

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
target_link_libraries(${PROJECT_NAME}_INTERFACE INTERFACE ${PROJECT_NAME} ${EXAGEOSTAT_DEP})
# Add linker options to the target
target_link_options(${PROJECT_NAME}_INTERFACE INTERFACE "SHELL:-Wl,--whole-archive $<TARGET_FILE:${PROJECT_NAME}> -Wl,--no-whole-archive")

# Install headers
install(TARGETS exageostatcpp
        DESTINATION lib/
        PUBLIC_HEADER DESTINATION include/
        )

# Add tests if enabled
if (${EXAGEOSTAT_BUILD_TESTS})
    message(STATUS "Building Tests")
    include(ImportCatch2)
    include(Catch)
    include(CTest)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp-tests)
    enable_testing()
endif ()


if (EXAGEOSTAT_BUILD_EXAMPLES)
    message(STATUS "Building Examples  is Enabled")
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/examples)
endif ()

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
message("\t export PKG_CONFIG_PATH=${EXAGEOSTAT_INSTALL_PREFIX}/CHAMELEON/lib/pkgconfig:${EXAGEOSTAT_INSTALL_PREFIX}/GSL/lib/pkgconfig:$PKG_CONFIG_PATH")
message("\t export PKG_CONFIG_PATH=${EXAGEOSTAT_INSTALL_PREFIX}/HCORE/lib/pkgconfig:${EXAGEOSTAT_INSTALL_PREFIX}/HICMA/lib/pkgconfig:$PKG_CONFIG_PATH")
message("\t export PKG_CONFIG_PATH=${EXAGEOSTAT_INSTALL_PREFIX}/HWLOC/lib/pkgconfig:${EXAGEOSTAT_INSTALL_PREFIX}/NLOPT/lib/pkgconfig:$PKG_CONFIG_PATH")
message("\t export PKG_CONFIG_PATH=${EXAGEOSTAT_INSTALL_PREFIX}/STARPU/lib/pkgconfig:${EXAGEOSTAT_INSTALL_PREFIX}/STARSH/lib/pkgconfig:$PKG_CONFIG_PATH")
message("\t ------------------------------------------------------------------------------------------------------------------------------- \n")
