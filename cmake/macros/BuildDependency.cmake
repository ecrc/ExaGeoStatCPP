# Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file BuildDependency.cmake
# @brief Fetches, builds, and installs a dependency.
# @version 1.0.0
# @author Mahmoud ElKarargy
# @author Sameh Abdulah
# @date 2023-03-12

# @param raw_name The name of the dependency.
# @param url The URL from which to fetch the dependency.
# @param tag The version or tag of the dependency to fetch.
# @param ${FLAGS} Additional flags to pass to the configure/make commands.
# @param ${ISCMAKE} A boolean flag indicating whether the dependency uses CMake as its build system.
# @param ${ISGIT} A boolean flag indicating whether the dependency is hosted on a git repository.

# This macro fetches the dependency using CMake's FetchContent module, and then builds and installs it.
# It also sets several environment variables (LD_LIBRARY_PATH, LIBRARY_PATH, CPATH, PKG_CONFIG_PATH,
# and ${capital_name}_DIR) and includes and links to the installation directory of the dependency.

# After building and installing the dependency, the macro installs the lib, include, and share directories  in the current directory.
macro(BuildDependency raw_name url tag ${FLAGS} ${ISCMAKE} ${ISGIT} ${AUTO_GEN})
    # Set the name of the dependency.
    string(TOLOWER ${raw_name} name)
    string(TOUPPER ${raw_name} capital_name)

    # Fetch the dependency, depending on whether it's a git repo or not.
    message(STATUS "Fetching ${name} ${tag} from ${url}")
    include(FetchContent)
    set(FETCHCONTENT_BASE_DIR ${PROJECT_SOURCE_DIR}/installdir/_deps/${capital_name}/)
    if (ISGIT)
        FetchContent_Declare(${name} GIT_REPOSITORY "${url}" GIT_TAG "${tag}")
    else()
        FetchContent_Declare(${name} URL "${url}")
    endif ()
    FetchContent_Populate(${name})

    # Set up build paths and create directory for build artifacts.
    set(${name}_srcpath ${PROJECT_SOURCE_DIR}/installdir/_deps/${capital_name}/${name}-src)
    set(${name}_binpath ${${name}_srcpath}/bin)
    set(${name}_installpath ${PROJECT_SOURCE_DIR}/installdir/_deps/${capital_name})
    file(MAKE_DIRECTORY ${${name}_binpath})

    # Configure subproject.
    if (ISCMAKE)
        execute_process(COMMAND ${CMAKE_COMMAND} ${FLAGS}
                ${${name}_srcpath}
                WORKING_DIRECTORY
                ${${name}_binpath})
    else()
        if (AUTO_GEN)
            execute_process(COMMAND ./autogen.sh
                    WORKING_DIRECTORY ${${name}_srcpath}
                    COMMAND_ERROR_IS_FATAL ANY)
        endif()
        execute_process(COMMAND ./configure ${FLAGS}
                WORKING_DIRECTORY ${${name}_srcpath}
                COMMAND_ERROR_IS_FATAL ANY)
    endif ()

    # Build and install subproject.
    include(ProcessorCount)
    ProcessorCount(N)
    if (ISCMAKE)
        execute_process(COMMAND make -j ${N}
                WORKING_DIRECTORY ${${name}_binpath}
                COMMAND_ERROR_IS_FATAL ANY)
        execute_process(COMMAND make install -j ${N}
                WORKING_DIRECTORY ${${name}_binpath}
                COMMAND_ERROR_IS_FATAL ANY)
    else()
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
    set(ENV{PKG_CONFIG_PATH} "${${name}_installpath}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
    set(${capital_name}_DIR "${${name}_installpath}")
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
