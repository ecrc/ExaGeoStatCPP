
# Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
# All rights reserved.
# ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

# @file ImportNlohmannJSON.cmake
# @brief Checks for the nlohmann library and includes it in the project if it is not already present.
# @version 1.1.0
# @author Mahmoud ElKarargy
# @date 2024-12-25

# Configurations for integrating the nlohmann
# 'name' is assigned "NlohmannJSON" to identify the nlohmann Library within this script.
set(name "nlohmann_json")
# 'tag' specifies the version tag "v3.11.2" for the nlohmann library, indicating the exact version to be used.
set(tag "v3.11.2")
# 'version' sets "3.11.2" as the version of the nlohmann library, ensuring compatibility with project requirements.
set(version "3.11.2")
# 'flag' is available for additional configuration options during build or installation, but remains empty here.
set(flag "")
# 'is_cmake' indicates whether nlohmann uses CMake for building.
set(is_cmake ON)
# 'is_git' denotes if nlohmann's source code is hosted in a Git repository.
set(is_git ON)
# 'auto_gen' signifies the need for autogen scripts in the build process. Here, it is set to OFF, indicating they are not needed.
set(auto_gen OFF)
# 'url' provides the download location for the nlohmann source code.
set(url "https://github.com/nlohmann/json.git")

# Include the 'ImportDependency' macro, responsible for managing the GSL library's import and setup process.
include(macros/ImportDependency)
# Execute the 'ImportDependency' macro with the previously established parameters to handle the detection, downloading, and integration of nlohmann.
ImportDependency(${name} ${tag} ${version} ${url} "${flag}" "" ${is_cmake} ${is_git} ${auto_gen})

# Add nlohmann_json to the project's list of linked libraries, making its functionality accessible within the project.
list(APPEND LIBS nlohmann_json::nlohmann_json)

# Output a message signaling the successful integration of the nlohmann library into the project.
message(STATUS "${name} done")
