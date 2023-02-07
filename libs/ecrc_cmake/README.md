ECRC CMake modules
====================

This project provides a collection of CMake modules that can be shared
among projects using CMake as build system.

For now it is mainly constituted of "Find" modules that help detecting
installed libraries on the system. These modules are located in

Get ecrc_cmake
---------------------

To use latest development states of ecrc_cmake, please clone the
master branch:

    git clone https://github.com/ecrc/ecrc_cmake.git

Documentation
---------------------

TODO

Installation
---------------------

To use ECRC modules you have to add the path to the modules in your
CMake project and include the EcrcInit module:

    # Define where are located module files on your system
    set(ECRC_CMAKE_MODULE_PATH "/where/is/ecrc_cmake" CACHE PATH "Path to ecrc_cmake sources")
    # Append this directory to the list of directories containing CMake modules
    list(APPEND CMAKE_MODULE_PATH "${ECRC_CMAKE_MODULE_PATH}/modules/" )
    # Include the init module
    include(EcrcInit)
    #

We recommend to use this project as a `git submodule` of your project.

Get involved!
---------------------

### Mailing list

TODO

### Contributions

https://github.com/ecrc/ecrc_cmake/blob/master/CONTRIBUTING.md

### Authors

The following people contributed to the development of ecrc_modules:
  * Cedric Castagnede
  * Mathieu Faverge, PI
  * Florent Pruvost, PI

If we forgot your name, please let us know that we can fix that mistake.

### Licence

https://github.com/ecrc/ecrc_cmake/blob/master/LICENCE.txt
