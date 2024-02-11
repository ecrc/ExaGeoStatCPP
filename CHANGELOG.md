# Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased] - Version 1.0.1 (YYYY-MM-DD)
### Added
- Implemented a new changelog.
- Introduced a benchmarking script.
- .gitignore file.
- More examples.
- Add tests for all src files.

### Fixed
- Resolved issues with MPI installation.
- Fixed the printing of configuration summaries with MPI.
- Automated the process of adding a new kernel.
- Improved packaging of software with CPack.
- Addressed installation issues.
- Fixed bivariate and trivariate kernels functionality.
- Corrected time-space kernel issues.

### Changed
- Updated the installation process for dependencies.
- Modified the calculation of P for all kernels.
- Adjusted CMake variables.
- Revised the process of finding BLASPP and Catch2 libraries.
- Updated doxygen documentation.
- Split the synthetic generator functions into BitHelper class and Locations generator class.
- Created a Bassel Function helper for kernels.
- Cleaned the code base for better readability.

### Removed
- Eliminated non-stationary kernel support.
- Removed Find OpenMP, LAPACKPP, and CuSOLVER.

## [1.0.0] - 2023-11-12
### Added
- Integrated all features present in [ExaGeoStat C version](https://github.com/ecrc/exageostat).
