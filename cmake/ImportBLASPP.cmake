# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file ImportBLASPP.cmake
# @brief This file searches for the BLAS++ library and includes it if not already included.
# @version 1.1.0
# @author Mahmoud ElKarargy
# @date 2024-02-04

# Include the ImportBLAS script, which presumably sets up basic BLAS dependencies.
include(ImportBLAS)

# Set up basic configuration variables for locating or installing BLAS++.
# `name` variable set to 'blaspp' to represent the BLAS++ library throughout the script.
set(name blaspp)
# Convert the name to uppercase and store it in `capital_name` for display purposes.
string(TOUPPER ${name} capital_name)
# Set the `tag` variable to specify the version of BLAS++ to be used or installed.
set(tag "v2023.01.00")
# Set the `url` variable to the location of the BLAS++ repository on GitHub.
set(url "https://github.com/icl-utk-edu/blaspp")
# Configure the directory where BLAS++ will be installed or found.
set(${name}_DIR "${CMAKE_INSTALL_PREFIX}/${capital_name}/${name}-build/")

# Display a header message indicating the start of the BLAS++ setup process.
message("")
message("---------------------------------------- ${capital_name}")
# Output the version of BLAS++ being checked for or installed.
message(STATUS "Checking for ${capital_name} with Version ${version}")

# Check if the BLAS++ target is already defined in the current CMake environment.
IF (NOT TARGET ${name})
    # Include the FindPkgConfig module to use pkg-config for finding installed libraries.
    include(FindPkgConfig)
    # Attempt to find the pkg-config utility quietly.
    find_package(PkgConfig QUIET)
    # Quietly search for the BLAS++ package, without specifying components or version.
    find_package(${name} ${version} QUIET COMPONENTS ${components})

    # If BLAS++ is found, output its location and libraries linked.
    if (${name}_FOUND)
        message("   Found ${capital_name}; ${${name}_DIR} ${${name}_LIBRARIES}")
    else ()
        # If BLAS++ is not found, indicate it and proceed to install it.
        message("   Can't find ${capital_name}, Installing it instead ..")
        # Include the FetchContent module to download external projects during the configure stage.
        include(FetchContent)
        # Set the base directory for downloading and building the external project.
        set(FETCHCONTENT_BASE_DIR ${CMAKE_INSTALL_PREFIX}/${capital_name})
        # Declare the external project (BLAS++) with its Git repository and specific tag.
        FetchContent_Declare(${name} GIT_REPOSITORY "${url}" GIT_TAG "${tag}")
        # Make the declared content available, effectively downloading and setting up BLAS++.
        FetchContent_MakeAvailable(${name})

    endif ()
else ()
    # If the BLAS++ target is already defined, indicate that it's already included.
    message(STATUS "${capital_name} already included")
endif ()

# Add BLAS++ to the list of libraries to be linked against by setting `LIBS`.
set(LIBS ${name} ${LIBS})
# Final status message indicating completion of the BLAS++ setup.
message(STATUS "${name} done")
