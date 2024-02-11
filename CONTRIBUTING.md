Contributing to ExaGeoStatCPP
--------------------------------------------------------------------------------

- [Synopsis](#synopsis)
- [Forking ExaGeoStatCPP](#forking-ExaGeoStatCPP)
- [C++ coding style](#c-coding-style)
- [Developing new features](#developing-new-features)
- [Developing bug fixes](#developing-bug-fixes)
- [Creating pull requests](#creating-pull-requests)
- [Tests](#tests)

Synopsis
--------------------------------------------------------------------------------

This guide is intended for developers seeking to contribute new functionalities or bug fixes to ExaGeoStatCPP.
It presumes a basic understanding of the `git` command-line tool and GitHub operations. 
The document will outline the criteria for a well-formed pull request and detail the requisite tests
that must be cleared for your contributions to be integrated into ExaGeoStatCPP.

C++ coding style
--------------------------------------------------------------------------------

Changes to ExaGeoStatCPP C/C++ code should conform to the
[Google C++ Style Guide](https://google.github.io/styleguide/cppguide.html) and
the following specific style details:

## Naming convention

### File Names

All file names should have a ```PascalCase```.

```c++
// Example
TemplateClass.hpp
TemplateClass.cpp
```

### Header Guards

All header guards should have a ```SCREAMING_SNAKE_CASE```.

```c++
// Example
#ifndef
TEMPLATE_ARCHITECTURE_TEMPLATE_CLASS_HPP
#define
TEMPLATE_ARCHITECTURE_TEMPLATE_CLASS_HPP

// Code

#endif //TEMPLATE_ARCHITECTURE_TEMPLATE_CLASS_HPP
```

### Namespaces

All namespaces should be all ```lowercase```.

```c++
// Example
namespace librarynamespace {

}
```

### Classes & Structs

All classes and structs should have a ```PascalCase```.

```c++
// Example

class TemplateClasss {

};
```

### Functions

All functions should have a ```PascalCase```.

```c++
// Example

void PublicFunction() {

}
```

### Variables

#### Public Variables

All public variables should have s ```PascalCase```.

```c++
// Example
int PublicVariable;
```

#### Members Variables

All members variables should start with ```m```followed by ```PascalCase```.

```c++
// Example
int mPrivateVariable;
```

#### Pointers Variables

All pointers variables should start with ```p``` followed by ```PascalCase```.

```c++
// Example
int* pPointerVariable;
```

#### Combinations

Private and pointer variables.

```c++
// Example
int* mpPrivatePointerVariable;
```

#### Argument Variables

All argument variables should start with ``a`` followed by ```PascalCase```. <b>If and only if the variable is an object
or pointer to object</b>

```c++
// Example

void Function(TemplateClass aTemplateClass) {

}

// For pointer case

void Function(TemplateClass apTemplateClass) {

}
```

#### Scope Variables

All scope variables should have a ```snake_case```.

```c++
// Example

int number_count;
```

## ExaGeoStat NameSpaces

- Default indentation is 4 spaces, and wrapped parameters have 4 spaces indent.

- In ExaGeoStatCPP, we have 4 namespaces, as follows:
    - `exageostat::api`: This namespace contains the high-level drivers for the ExaGeoStatCPP functionalities that are
      provided to library users. These functions help users interact with the ExaGeoStatCPP framework and perform
      various statistical operations.
    - `exageostat::common`: This namespace contains all ExaGeoStatCPP common functionalities that might be used across
      the different modules of the ExaGeoStatCPP framework.
    - `exageostat::configurations`: This namespace contains all ExaGeoStatCPP configurations arguments and parsers.
      These functions are used to parse and set the configuration parameters for the ExaGeoStatCPP framework.
    - `exageostat::data-generators`: This namespace is used to generate data sets.
    - `exageostat::data-units`: This namespace is used for all ExaGeoStatCPP base data structures that the user should
      utilize and interact with. These data units are used to represent the data and perform operations on it.
    - `exageostat::helpers`: This namespace contains helper functions that can be used across the different modules of
      the ExaGeoStatCPP framework.
    - `exageostat::kernels`: These functions provide low-level implementations of the supported kernels offered by the
      ExaGeoStatCPP framework.
    - `exageostat::linear-algebra-solvers`: This namespace is used for all ExaGeoStatCPP integrated linear algebra
      solvers libraries.
    - `exageostat::operators`: This namespace contains various operators used by the ExaGeoStatCPP framework. These
      operators are used to perform various mathematical operations on the data sets.

Use `clang-format` to check your C/C++ changes.

To install on Ubuntu 16+, do:

    $ apt-get install -y clang-format

You can check the format of a C/C++ file with the following:

    $ clang-format <path/to/file.cc> --style=google > path/to/file_marked.cc
    $ diff <path/to/file.cc> path/to/file_marked.cc

Forking ExaGeoStatCPP
--------------------------------------------------------------------------------

If you're not affiliated with KAUST, you won't possess the rights to push new branches
directly to the repository. Your initial step should be to create a fork of the ExaGeoStatCPP repository. 
This action will generate a version of the repository under your ownership, enabling you to upload your
changes to GitHub and initiate pull requests.

Developing new features
--------------------------------------------------------------------------------

When developing new features, you should base your work on the main branch.
To begin implementing new functionality, first, make sure your local main branch
is synchronized with the latest updates:

    $ git checkout main
    $ git pull origin main

You can now create a new branch to develop your feature on:

    $ git checkout -b feature/<name-of-feature>

Begin the development of your feature on this branch, ensuring to include tests that
validate your new code. If you're introducing new methods or classes, remember to
provide Doxygen documentation for these additions.

After completing your feature and confirming that your tests succeed, you may push your
branch to GitHub and initiate a pull request.

Developing bug fixes
--------------------------------------------------------------------------------

Initially, verify whether the modification you're aiming to implement has already been addressed
in the main branch. If it has, we recommend switching to using the main branch or temporarily
integrating the fix into your current version of ExaGeoStatCPP.

If the bug remains unresolved, ensure that you have the most recent version of the main branch:

    $ git checkout main
    $ git pull origin main

Then, create a new branch for your bugfix:

    $ git checkout -b bugfix/<name-of-bug>


Begin by adding a test that replicates the bug you've discovered. Proceed with developing
your bug fix, as usual, incorporates tests that verify your changes have indeed fixed the issue.

Upon completion, push your branch to GitHub and proceed to create a pull request.

Adding new kernel
--------------------------------------------------------------------------------
To add a new kernel, base your work on the main branch. Before you start developing the new feature,
make sure your local main branch is current with the latest updates:

    $ git checkout main
    $ git pull origin main

You can now create a new branch to develop your feature on:

    $ git checkout -b kernel/<name-of-kernel>

Proceed to develop your feature on this branch and add tests that will exercise
your new code. If you are creating new methods or classes, please add Doxygen
documentation.

To add a new kernel, you need to follow these steps.

1. Place your kernel header file in the inst/include/kernels/concrete directory. The file name should match the kernel's name. For instance, if your header file is named UnivariateMaternStationary.hpp, it can be invoked using either univariate_matern_stationary or UnivariateMaternStationary. The naming linkage is handled automatically, so there's no additional setup required on your part.

2. Derive from the base class located in Kernel.hpp and implement the necessary functions.

3. Ensure your kernel includes all the requisite functions that adhere to the established naming conventions found in other kernels. This will allow for proper support and integration of your new kernel.

After finalizing your feature and confirming that your tests run successfully, 
you are ready to push your branch to GitHub and submit a pull request.

Creating pull requests
--------------------------------------------------------------------------------

You can create a new pull request
[here](https://github.com/ecrc/ExaGeoStatCPP-dev/pulls).
Ensure that your pull request base is the `main` branch of ExaGeoStatCPP.

Add a descriptive title explaining the bug you fixed or the feature you have
added, and put a longer description of the changes you have made in the comment
box.

After you have submitted your pull request, it will undergo automated testing
and also receive a review by members of the HiCMA team. If your branch successfully
clears both the automated tests and the peer reviews, it will be approved for merging
into ExaGeoStatCPP.

Tests
--------------------------------------------------------------------------------

ExaGeoStatCPP employs Jenkins for continuous integration testing. Each pull request 
automatically triggers our suite of tests, and a prerequisite for merging your pull 
request is the successful passage of all these tests. When you're working on a bug 
fix or introducing a new feature, it is crucial to include tests that validate the 
integrity of your code. Given that ExaGeoStatCPP operates across an array of systems and 
configurations, introducing new tests is instrumental in guaranteeing that every 
feature functions correctly in diverse settings.

In ExaGeoStatCPP, the tests are organized within the tests directory, with each set 
of tests categorically divided by component.
