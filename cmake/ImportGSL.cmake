
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
        message("   Found GSL: ${GSL_INCLUDE_DIRS}")
    else ()

        # If the GSL library is not found, install it and add it to the project's libraries.
        message("   Can't find GSL, Installing it instead ..")
        set(FLAGS "")
        set(ISCMAKE OFF)
        set(ISGIT OFF)
        set(AUTO_GEN OFF)
        BuildDependency(GSL "https://ftp.gnu.org/gnu/gsl/gsl-2.7.1.tar.gz" "v2.7.1" ${FLAGS} ${ISCMAKE} ${ISGIT} ${AUTO_GEN})
        find_package(GSL REQUIRED)
    endif ()
else ()
    message("   GSL already included")
endif ()

# Add the GSL library to the project's list of libraries.
list(APPEND LIBS gsl)
include_directories(${GSL_INCLUDE_DIRS})
link_directories(${GSL_LIBRARY_DIRS})
message("${GSL_LIBRARIES}")

message(STATUS "GSL done")
