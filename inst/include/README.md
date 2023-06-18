# Include Subdirectory
This contains all the header files for the project.

## File structure
- `api`: Directory contains the high-level drivers for the ExaGeoStat-cpp functionalities that are provided to library users. These functions help users interact with the ExaGeoStat-cpp framework and perform various statistical operations.
- `common`: Directory contains all ExaGeoStat-cpp common and helper functionalities that might provide useful facilities in the examples and testing. These functions provide common functionality that can be used across the different modules of the ExaGeoStat-cpp framework.
- `configurations`: Directory contains all ExaGeoStat-cpp configurations arguments and parsers. These functions are used to parse and set the configuration parameters for the ExaGeoStat-cpp framework.
- `data-generators`: Directory contains is used for ExaGeoStat-cpp implementations of the various options for data-generation module. These functions generate synthetic data sets that can be used for testing and demonstration purposes.
- `data-units`: Directory contains is used for all ExaGeoStat-cpp base data structures that the user should utilize and interact with, including the base tiles. These data units are used to represent the data and perform operations on it.
- `helpers`: Directory contains helper functions used by other modules of the ExaGeoStat-cpp framework. These functions provide common functionality that can be used across the different modules of the ExaGeoStat-cpp framework.
- `kernels`: Directory contains is used for target backend implementations of the various kernels used internally to fulfill the targeted operations. These functions provide low-level implementations of the operations performed by the ExaGeoStat-cpp framework.
- `linear-algebra-solvers`: Directory contains is used for all ExaGeoStat-cpp integrated linear algebra solvers libraries. These solvers are used to solve the linear algebra problems that arise during the execution of the ExaGeoStat-cpp framework.
- `operators`: Directory contains various operators used by the ExaGeoStat-cpp framework. These operators are used to perform various mathematical operations on the data sets.