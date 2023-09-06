
# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file ImportGSL.cmake
# @brief Checks for the GSL library and includes it in the project if it is not already present.
# @version 1.0.0
# @author Mahmoud ElKarargy
# @author Sameh Abdulah
# @date 2023-03-16

message("")
message("---------------------------------------- GSL")
message(STATUS "Checking for GSL")

include(macros/BuildDependency)

# Check if the GSL library is already included in the project.
if (NOT TARGET GSL_FOUND)

    # If not, attempt to find it with PkgConfig and CMake's find_package function.
    include(FindPkgConfig)
    find_package(PkgConfig QUIET)
    find_package(GSL QUIET)

    # If the GSL library is found, add it to the project's libraries.
    if (GSL_FOUND)
        message("   Found GSL: ${GSL_LIBDIR}")
    else()

        # If the GSL library is not found, install it and add it to the project's libraries.
        message("   Can't find GSL, Installing it instead ..")
        set(FLAGS --prefix=${PROJECT_SOURCE_DIR}/installdir/_deps/GSL/)
        set(ISCMAKE OFF)
        set(ISGIT OFF)
        set(AUTO_GEN OFF)
        BuildDependency(GSL "https://ftp.gnu.org/gnu/gsl/gsl-2.6.tar.gz" "v2.6" ${FLAGS} ${ISCMAKE} ${ISGIT} ${AUTO_GEN})
        set(FLAGS "")
        find_package(GSL REQUIRED)
    endif()
else()
    message("   GSL already included")
endif()

# Add the GSL library to the project's list of libraries.
list(APPEND LIBS gsl)
set(ENV{CPATH} "${PROJECT_SOURCE_DIR}/installdir/_deps/GSL/include:$ENV{CPATH}")
include_directories(${PROJECT_SOURCE_DIR}/installdir/_deps/GSL/include)
link_directories(${PROJECT_SOURCE_DIR}/installdir/_deps/GSL/lib)

message(STATUS "GSL done")