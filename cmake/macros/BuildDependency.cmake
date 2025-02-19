# Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file BuildDependency.cmake
# @brief Fetches, builds, and installs a dependency.
# @version 1.1.0
# @author Mahmoud ElKarargy
# @author Amr Nasr
# @date 2024-02-04

# After building and installing the dependency, the macro installs the lib, include, and share directories
# in the current directory.

# BuildDependency Macro:
# This macro is designed to fetch, configure, build, and install a dependency.
# It takes the following parameters:
# - raw_name: The name of the dependency.
# - url: The URL of the repository or source tarball.
# - tag: The version or tag of the dependency to fetch.
# - flags: Additional flags to pass to the configure/make commands.
# - is_using_cmake: A boolean flag indicating whether the dependency uses CMake as its build system.
# - is_using_git: A boolean flag indicating whether the dependency is hosted on a git repository.
# - auto_generation: A boolean flag indicating whether to use autogen scripts or not.

# The macro fetches the dependency using CMake's FetchContent module, depending on whether it's a git repo or not.
# It sets up build paths and creates a directory for build artifacts. The subproject is then configured using
# CMake or autotools, and finally, it's built and installed. Environment variables are set, and the dependency's
# lib, include, and share directories are installed in the current directory.

# Define the macro BuildDependency with parameters for handling various aspects of dependency management.
macro(BuildDependency raw_name url tag flags is_using_cmake is_using_git auto_generation)

    # Convert the raw dependency name to lowercase and uppercase for different uses and set them as 'name' and 'capital_name'.
    string(TOLOWER ${raw_name} name)
    string(TOUPPER ${raw_name} capital_name)

    # Log the start of the fetch process for the dependency, including its name, tag, and source URL.
    message(STATUS "Fetching ${name} ${tag} from ${url}")
    # Include the CMake module for downloading and updating content during the configure step.
    include(FetchContent)
    # Set the base directory for fetched content to a directory within the install prefix, named after the dependency.
    set(FETCHCONTENT_BASE_DIR ${CMAKE_INSTALL_PREFIX}/${capital_name})

    # Check if the dependency is hosted in a git repository and declare it accordingly with FetchContent, using git-specific options.
    if (${is_using_git})
        FetchContent_Declare(${name}
                GIT_REPOSITORY "${url}"
                GIT_TAG "${tag}"
                GIT_SHALLOW TRUE  # For a shallow clone, fetching only the history needed for the specified tag
                GIT_PROGRESS TRUE  # Show progress during the clone
                )
    else ()
        # If not using git, declare the dependency for FetchContent using a direct URL (e.g., for a tarball).
        FetchContent_Declare(${name} URL "${url}")
    endif ()

    # Make the content available, effectively downloading it if necessary.
    FetchContent_Populate(${name})

    # Set variables for the source path, binary (build) path, and installation path of the dependency.
    set(${name}_srcpath ${CMAKE_INSTALL_PREFIX}/${capital_name}/${name}-src)
    set(${name}_binpath ${${name}_srcpath}/bin)
    set(${name}_installpath ${CMAKE_INSTALL_PREFIX}/${capital_name})
    # Ensure the binary path directory exists.
    file(MAKE_DIRECTORY ${${name}_binpath})

    # Configure the project. If using CMake, run cmake command with specified flags and install prefix within the binary path.
    if (${is_using_cmake})
        execute_process(COMMAND ${CMAKE_COMMAND} -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}/${capital_name} -DCMAKE_C_FLAGS=-fPIC ${flags}
                ${${name}_srcpath}
                WORKING_DIRECTORY ${${name}_binpath})
    else ()
        # For non-CMake projects, run autogen.sh if auto_generation is true, then configure the project with specified flags.
        if (${auto_generation})
            execute_process(COMMAND ./autogen.sh
                    WORKING_DIRECTORY ${${name}_srcpath}
                    COMMAND_ERROR_IS_FATAL ANY)  # Halt on error
        endif ()
        execute_process(COMMAND ./configure --prefix=${CMAKE_INSTALL_PREFIX}/${capital_name} ${flags}
                WORKING_DIRECTORY ${${name}_srcpath}
                COMMAND_ERROR_IS_FATAL ANY)  # Halt on error
    endif ()

    # Include the ProcessorCount module to determine the number of CPUs for parallel build and install commands.
    include(ProcessorCount)
    ProcessorCount(N)
    # Subtract 5 from N, ensuring it doesn't go below 0
    math(EXPR N "${N} - 5")
    if (N LESS 0)
        set(N 1)
    endif()
    # Build the project using make, with parallel jobs based on processor count. This applies to both CMake and non-CMake projects.
    if (${is_using_cmake})
        execute_process(COMMAND make -j ${N}
                WORKING_DIRECTORY ${${name}_binpath}
                COMMAND_ERROR_IS_FATAL ANY)  # Halt on error
        # Install the built project, also with parallel jobs.
        execute_process(COMMAND make install -j ${N}
                WORKING_DIRECTORY ${${name}_binpath}
                COMMAND_ERROR_IS_FATAL ANY)  # Halt on error
    else ()
        execute_process(COMMAND make -j ${N}
                WORKING_DIRECTORY ${${name}_srcpath}
                COMMAND_ERROR_IS_FATAL ANY)  # Halt on error
        execute_process(COMMAND make install -j ${N}
                WORKING_DIRECTORY ${${name}_srcpath}
                COMMAND_ERROR_IS_FATAL ANY)  # Halt on error
    endif ()

    # Set environment variables for dynamic and static linking as well as include paths, pointing to the dependency's installation directory.
    set(ENV{LD_LIBRARY_PATH} "${${name}_installpath}/lib:${${name}_installpath}/lib64:$ENV{LD_LIBRARY_PATH}")
    set(ENV{LIBRARY_PATH} "${${name}_installpath}/lib:${${name}_installpath}/lib64:$ENV{LIBRARY_PATH}")
    set(ENV{CPATH} "${${name}_installpath}/include:$ENV{CPATH}")
    set(ENV{PKG_CONFIG_PATH} "${${name}_installpath}/lib/pkgconfig:${${name}_installpath}/lib64/pkgconfig:$ENV{PKG_CONFIG_PATH}")

    # Include and link to the installation directory of the dependency
    include_directories(${${name}_installpath}/include)
    link_directories(${${name}_installpath}/lib)

    # Install the dependency's lib, include, and share directories in the current directory.
    install(
            DIRECTORY
            "${${name}_installpath}/lib"
            DESTINATION
            .
    )
    install(
            DIRECTORY
            "${${name}_installpath}/include"
            DESTINATION
            .
    )
    install(
            DIRECTORY
            "${${name}_installpath}/share"
            DESTINATION
            .
    )
endmacro()