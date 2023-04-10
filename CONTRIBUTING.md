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

All pointers variables should start with ```p```followed by ```PascalCase```.

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

All argument variables should start with ``a``followed by ```PascalCase```. <b>If and only if the variable is an object
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

- In ExaGeoStat-CPP, we have 4 namespaces, as follows:
    - `exageostat::common`: Namespace used for all ExaGeoStat-CPP common and helper functionalities that might provide useful facilities in the examples and testing.
    - `exageostat::configurations`: Namespace used for ExaGeoStat-CPP configurations arguments and parsers.
    - `exageostat::kernels`: Namespace used for target backend implementations of the various kernels used internally to fulfill the targeted operations.
    - `exageostat::data-generators`: Namespace used for ExaGeoStat-CPP implementations of the various options for data-generation module.
    - `exageostat::data-units`: Namespace used for all ExaGeoStat-CPP base data structures that the user should utilize and interact with, including the base tiles.
    - `exageostat::linear-algebra-solvers`: Namespace used for all ExaGeoStat-CPP integrated linear algebra solvers libraries.

Use `clang-format` to check your C/C++ changes.

To install on Ubuntu 16+, do:

    $ apt-get install -y clang-format

You can check the format of a C/C++ file with the following:

    $ clang-format <path/to/file.cc> --style=google > path/to/file_marked.cc
    $ diff <path/to/file.cc> path/to/file_marked.cc

Forking ExaGeoStat-CPP
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

    $ git checkout -b FT-<name-of-feature>-<FIRST-THREE-LETTERS-OF-YOUR-NAME>

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

    $ git checkout -b FIX-<name-of-bug>-<FIRST-THREE-LETTERS-OF-YOUR-NAME>

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
