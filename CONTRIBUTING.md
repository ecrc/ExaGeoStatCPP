Contributing to ExaGeoStat-CPP
--------------------------------------------------------------------------------

- [Synopsis](#synopsis)
- [Forking ExaGeoStat-CPP](#forking-ExaGeoStat-CPP)
- [C++ coding style](#c-coding-style)
- [Developing new features](#developing-new-features)
- [Developing bug fixes](#developing-bug-fixes)
- [Creating pull requests](#creating-pull-requests)
- [Tests](#tests)

Synopsis
--------------------------------------------------------------------------------

This document is intended for developers who want to add new features or
bug fixes to ExaGeoStat-CPP. It assumes you have some familiarity with git and GitHub. It
will discuss what a good pull request looks like, and the tests that your
pull request must pass before it can be merged into ExaGeoStat-CPP.

C++ coding style
--------------------------------------------------------------------------------

Changes to ExaGeoStat-CPP C/C++ code should conform to the
[Google C++ Style Guide](https://google.github.io/styleguide/cppguide.html) and
the following HiCMA, HCore specific style details:

- Naming convention
    - Global functions are snake_case.
    - Member functions are UpperCamelCase.
    - Classes, enums, enum values, structs are UpperCamelCase.
    - Variables are snake_case.
    - Constants are snake_case.
    - Namespaces are lowercase.
    - Member private variables have an m suffix; Note that using leading
      underscores is illegal (they are reserved by C/C++).
- Default indentation is 4 spaces, and wrapped parameters have 4 spaces indent.
- In ExaGeoStat-CPP, we have 4 namespaces, as follows:
    - `exageostat::dataunits`: Namespace used for ExaGeoStat-CPP main memory unit used to contain some contiguous sub-matrix elements.
    - `exageostat::kernels`: Namespace used for target backend implementations of the various kernels used internally to fulfill the targeted operations.
    - `exageostat::cudakernels`: Namespace used for CUDA backend implementations of the various kernels used internally to fulfill the targeted operations.
    - `exageostat::operators`: Namespace used for all ExaGeoStat-CPP base data structures that the user should utilize and interact with, including the base tiles.
    - `exageostat::helpers`: Namespace used for all ExaGeoStat-CPP helper functionalities that might provide useful facilities in the examples and testing, but not critical for ExaGeoStat-CPP core operations.

Use `clang-format` to check your C/C++ changes.

To install on Ubuntu 16+, do:

    $ apt-get install -y clang-format

You can check the format of a C/C++ file with the following:

    $ clang-format <path/to/file.cc> --style=google > path/to/file_marked.cc
    $ diff <path/to/file.cc> path/to/file_marked.cc

Forking HExaGeoStat-CPP
--------------------------------------------------------------------------------

If you aren't a HiCMA or HCore developer at KAUST, then you won't have permission to push
new branches to the repository. First, you should create a fork. This will
create a copy of the ExaGeoStat repository that you own, and will ensure you can push
your changes up to GitHub and create pull requests.

Developing new features
--------------------------------------------------------------------------------

New features should be based on the `main` branch. When you want to create a
new feature, first ensure you have an up-to-date copy of the `main` branch:

    $ git checkout main
    $ git pull origin main

You can now create a new branch to develop your feature on:

    $ git checkout -b feature/<name-of-feature>

Proceed to develop your feature on this branch, and add tests that will exercise
your new code. If you are creating new methods or classes, please add Doxygen
documentation.

Once your feature is complete and your tests are passing, you can push your
branch to GitHub and create a pull request.

Developing bug fixes
--------------------------------------------------------------------------------

First, check if the change you want to make has been fixed in `main`. If so,
we suggest you either start using the `main` branch, or temporarily apply the
fix to whichever version of ExaGeoStat-CPP you are using.

Assuming there is an unsolved bug, first make sure you have an up-to-date copy
of the `main` branch:

    $ git checkout main
    $ git pull origin main

Then create a new branch for your bugfix:

    $ git checkout -b bugfix/<name-of-bug>

First, add a test that reproduces the bug you have found. Then develop your
bugfix as normal, and add tests to check your changes actually fix the bug.

Once you are finished, you can push your branch to GitHub, then create a pull
request.

Creating pull requests
--------------------------------------------------------------------------------

You can create a new pull request
[here](https://github.com/ecrc/ExaGeoStat-CPP-dev/pulls).
Ensure that your pull request base is the `main` branch of ExaGeoStat-CPP.

Add a descriptive title explaining the bug you fixed or the feature you have
added, and put a longer description of the changes you have made in the comment
box.

Once your pull request has been created, it will be run through our automated
tests and also be reviewed by HiCMA team members. Providing the branch passes
both the tests and reviews, it will be merged into ExaGeoStat-CPP.

Tests
--------------------------------------------------------------------------------

ExaGeoStat-CPP uses Jenkins for continuous integration tests. Our tests are automatically
run against every new pull request, and passing all tests is a requirement for
merging your pull request. If you are developing a bugfix or a new feature,
please add a test that checks the correctness of your new code. ExaGeoStat-CPP is used on
a wide variety of systems with a number of configurations, and adding new tests
helps ensure that all features work as expected across these environments.

ExaGeoStat-CPP's tests are all in the `tests` directory and are split up by component.
