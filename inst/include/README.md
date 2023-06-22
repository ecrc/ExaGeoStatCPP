# Include Subdirectory
This contains all the header files for the project.

## File structure
- `api`: Directory contains the high-level drivers for the ExaGeoStat-cpp functionalities that are provided to library users. These functions help users interact with the ExaGeoStat-cpp framework and perform various statistical operations.
- `common`: Directory contains all ExaGeoStat-cpp common functionalities that might be used across the different modules of the ExaGeoStat-cpp framework.
- `configurations`: Directory contains all ExaGeoStat-cpp configurations arguments and parsers. These functions are used to parse and set the configuration parameters for the ExaGeoStat-cpp framework.
- `data-generators`: Directory contained the needed methods to generate data sets.
- `data-units`: Directory is used for all ExaGeoStat-cpp base data structures that the user should utilize and interact with. These data units are used to represent the data and perform operations on it.
- `helpers`: Directory contains helper functions that can be used across the different modules of the ExaGeoStat-cpp framework.
- `kernels`: Directory provide low-level implementations of the supported kernels offered by the ExaGeoStat-cpp framework.
- `linear-algebra-solvers`: Directory is used for all ExaGeoStat-cpp integrated linear algebra solvers libraries.
- `operators`: Directory contains various operators used by the ExaGeoStat-cpp framework. These operators are used to perform various mathematical operations on the data sets.
