# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file BuildDependency.cmake
# @brief Fetches, builds, and installs a dependency.
# @version 1.0.1
# @author Mahmoud ElKarargy
# @date 2023-03-12

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

macro(BuildDependency raw_name url tag flags is_using_cmake is_using_git auto_generation)

    # Set the name of the dependency.
    string(TOLOWER ${raw_name} name)
    string(TOUPPER ${raw_name} capital_name)

    # Fetch the dependency, depending on whether it's a git repo or not.
    message(STATUS "Fetching ${name} ${tag} from ${url}")
    include(FetchContent)
    set(FETCHCONTENT_BASE_DIR ${CMAKE_INSTALL_PREFIX}/${capital_name})

    if (${is_using_git})
        FetchContent_Declare(${name}
                GIT_REPOSITORY "${url}"
                GIT_TAG "${tag}"
                GIT_SHALLOW TRUE
                GIT_PROGRESS TRUE
                )
    else ()
        FetchContent_Declare(${name} URL "${url}")
    endif ()

    FetchContent_Populate(${name})

    # Set up build paths and create a directory for build artifacts.
    set(${name}_srcpath ${CMAKE_INSTALL_PREFIX}/${capital_name}/${name}-src)
    set(${name}_binpath ${${name}_srcpath}/bin)
    set(${name}_installpath ${CMAKE_INSTALL_PREFIX}/${capital_name})
    file(MAKE_DIRECTORY ${${name}_binpath})

    # Configure subproject.
    if (${is_using_cmake})
        execute_process(COMMAND ${CMAKE_COMMAND} -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}/${capital_name} -DCMAKE_C_FLAGS=-fPIC ${flags}
                ${${name}_srcpath}
                WORKING_DIRECTORY ${${name}_binpath})
    else ()
        if (${auto_generation})
            execute_process(COMMAND ./autogen.sh
                    WORKING_DIRECTORY ${${name}_srcpath}
                    COMMAND_ERROR_IS_FATAL ANY)
        endif ()
        execute_process(COMMAND ./configure --prefix=${CMAKE_INSTALL_PREFIX}/${capital_name} ${flags}
                WORKING_DIRECTORY ${${name}_srcpath}
                COMMAND_ERROR_IS_FATAL ANY)
    endif ()

    # Build and install subproject.
    include(ProcessorCount)
    ProcessorCount(N)
    if (${is_using_cmake})
        execute_process(COMMAND make -j ${N}
                WORKING_DIRECTORY ${${name}_binpath}
                COMMAND_ERROR_IS_FATAL ANY)
        execute_process(COMMAND make install -j ${N}
                WORKING_DIRECTORY ${${name}_binpath}
                COMMAND_ERROR_IS_FATAL ANY)
    else ()
        execute_process(COMMAND make -j ${N}
                WORKING_DIRECTORY ${${name}_srcpath}
                COMMAND_ERROR_IS_FATAL ANY)
        execute_process(COMMAND make install -j ${N}
                WORKING_DIRECTORY ${${name}_srcpath}
                COMMAND_ERROR_IS_FATAL ANY)
    endif ()

    # Set environment variables and include/link to the installation directory of the dependency.
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