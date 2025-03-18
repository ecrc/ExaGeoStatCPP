
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file Definitions.hpp
 * @version 1.1.0
 * @brief This file contains common definitions used in ExaGeoStat software package.
 * @details These definitions include enums for dimension, computation, precision, and floating point arithmetic;
 * A macro for instantiating template classes with supported types; and a set of available kernels.
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-03-21
**/

#ifndef EXAGEOSTATCPP_DEFINITIONS_HPP
#define EXAGEOSTATCPP_DEFINITIONS_HPP

// This is a fix to avoid the problem in HiCMA which set the min definition with a conflict implementation of chrono library.
#ifdef USE_HICMA
#ifdef min
#undef min
#endif
#endif

#include <set>
#include <filesystem>

/**
 * @def EXAGEOSTAT_INSTANTIATE_CLASS
 * @brief Macro definition to instantiate the EXAGEOSTAT template classes with supported types.
 *
**/
#define EXAGEOSTAT_INSTANTIATE_CLASS(TEMPLATE_CLASS)   template class TEMPLATE_CLASS<float>;  \
                                                    template class TEMPLATE_CLASS<double>;

// Variables sizes.
#define SIZE_OF_FLOAT 4
#define SIZE_OF_DOUBLE 8

/**
 * Pi value.
 */
#define PI (3.141592653589793)

/**
 * Earth Radius value.
 */
#define EARTH_RADIUS 6371.0

/**
 * Q Norm value.
 */
#define Q_NORM 1.959964

/**
 * Logging Path Definition
 */
#define LOG_PATH PROJECT_SOURCE_DIR "/synthetic_ds/"

namespace exageostat::common {

    /**
     * @enum VerbosityLevel
     * @brief Enum denoting the run mode
     *
     */
    enum Verbose {
        QUIET_MODE = 0,
        STANDARD_MODE = 1,
        DETAILED_MODE = 2
    };

    /**
     * @enum Dimension
     * @brief Enum denoting the dimension of generated data.
     *
     */
    enum Dimension {
        Dimension2D = 0,
        Dimension3D = 1,
        DimensionST = 2,
    };

    /**
     * @enum Side
     * @brief Enum denoting the side on which the matrix appears in an equation.
     *
     */
    enum Side {
        EXAGEOSTAT_LEFT = 141,
        EXAGEOSTAT_RIGHT = 142
    };

    /**
     * @enum Trans
     * @brief Enum denoting whether or not to transpose a matrix.
     *
     */
    enum Trans {
        EXAGEOSTAT_NO_TRANS = 111,
        EXAGEOSTAT_TRANS = 112,
        EXAGEOSTAT_CONJ_TRANS = 113
    };

    /**
     * @enum Diag
     * @brief Enum denoting whether the diagonal is unitary.
     *
     */
    enum Diag {
        EXAGEOSTAT_NON_UNIT = 131,
        EXAGEOSTAT_UNIT = 132
    };

    /**
     * @enum Distance metric
     * @brief Enum denoting distance metric type.
     *
     */
    enum DistanceMetric {
        EUCLIDEAN_DISTANCE = 0,
        GREAT_CIRCLE_DISTANCE = 1
    };

    /**
     * @enum Descriptor Type
     * @brief Enum denoting the Descriptor Type.
     *
     */
    enum DescriptorType {
        CHAMELEON_DESCRIPTOR = 0,
        HICMA_DESCRIPTOR = 1,
        PARSEC_DESCRIPTOR = 2
    };

    /**
     * @enum Data source Type
     * @brief Enum denoting the data source Type.
     *
     */
    enum DataSourceType {
        SYNTHETIC = 0,
        CSV_FILE = 1,
        PARSEC_FILE = 2
    };

    /**
     * @enum Descriptor Name
     * @brief Enum denoting all Descriptors Names.
     *
     */
    enum DescriptorName : int {
        DESCRIPTOR_C = 0,
        DESCRIPTOR_Z = 1,
        DESCRIPTOR_Z_COPY = 2,
        DESCRIPTOR_PRODUCT = 3,
        DESCRIPTOR_DETERMINANT = 4,
        DESCRIPTOR_CD = 5,
        DESCRIPTOR_CUV = 6,
        DESCRIPTOR_CRK = 7,
        DESCRIPTOR_Z_OBSERVATIONS = 8,
        DESCRIPTOR_Z_Actual = 9,
        DESCRIPTOR_Z_MISS = 10,
        DESCRIPTOR_MSPE = 11,
        DESCRIPTOR_Z_1 = 12,
        DESCRIPTOR_Z_2 = 13,
        DESCRIPTOR_Z_3 = 14,
        DESCRIPTOR_PRODUCT_1 = 15,
        DESCRIPTOR_PRODUCT_2 = 16,
        DESCRIPTOR_PRODUCT_3 = 17,
        DESCRIPTOR_C11 = 18,
        DESCRIPTOR_C12 = 19,
        DESCRIPTOR_C22 = 20,
        DESCRIPTOR_C12D = 21,
        DESCRIPTOR_C12UV = 22,
        DESCRIPTOR_C12RK = 23,
        DESCRIPTOR_C22D = 24,
        DESCRIPTOR_C22UV = 25,
        DESCRIPTOR_C22RK = 26,
        DESCRIPTOR_MSPE_1 = 27,
        DESCRIPTOR_MSPE_2 = 28,
        DESCRIPTOR_k_T = 29,
        DESCRIPTOR_k_A = 30,
        DESCRIPTOR_k_A_TMP = 31,
        DESCRIPTOR_k_T_TMP = 32,
        DESCRIPTOR_K_T = 33,
        DESCRIPTOR_K_T_TMP = 34,
        DESCRIPTOR_K_A = 35,
        DESCRIPTOR_EXPR_1 = 36,
        DESCRIPTOR_EXPR_2 = 37,
        DESCRIPTOR_EXPR_3 = 38,
        DESCRIPTOR_EXPR_4 = 39,
        DESCRIPTOR_MLOE = 40,
        DESCRIPTOR_MMOM = 41,
        DESCRIPTOR_MLOE_MMOM = 42,
        DESCRIPTOR_ALPHA = 43,
        DESCRIPTOR_TRUTH_ALPHA = 44,
        DESCRIPTOR_TIMATED_ALPHA = 45,
        DESCRIPTOR_CK = 46,
        DESCRIPTOR_CJ = 47,
        DESCRIPTOR_C_TRACE = 48,
        DESCRIPTOR_C_DIAG = 49,
        DESCRIPTOR_A = 50,
        DESCRIPTOR_RESULTS = 51,
        DESCRIPTOR_SUM = 52,
        DESCRIPTOR_R = 53,
        DESCRIPTOR_R_COPY = 54,
        DESCRIPTOR_F_DATA = 55,
        DESCRIPTOR_ET1 = 56,
        DESCRIPTOR_ET2 = 57,
        DESCRIPTOR_EP = 58,
        DESCRIPTOR_SLMN = 59,
        DESCRIPTOR_IE = 60,
        DESCRIPTOR_IO = 61,
        DESCRIPTOR_P = 62,
        DESCRIPTOR_D = 63,
        DESCRIPTOR_FLMERA = 64,
        DESCRIPTOR_ZLM = 65,
        DESCRIPTOR_SC = 66,
        DESCRIPTOR_F_SPATIAL = 67,
        DESCRIPTOR_FLM = 68,
        DESCRIPTOR_FLMT = 69
    };

    /**
     * @enum Tile Storage
     * @brief Matrix tile storage
     *
     */
    typedef enum TileStorage {
        EXAGOSTAT_CM = 101,
        EXAGOSTAT_RM = 102,
        EXAGOSTAT_CCRB = 103,
        EXAGOSTAT_CRRB = 104,
        EXAGOSTAT_RCRB = 105,
        EXAGOSTAT_RRRB = 106,
    } ExaGeoStatTileStorage;

    /**
     * @enum Computation
     * @brief Enum denoting the types of computations that can be requested,
     * to use the required Linear Algebra solver library.
     *
     */
    enum Computation {
        EXACT_DENSE = 0,
        DIAGONAL_APPROX = 1,
        TILE_LOW_RANK = 2,
    };

    /**
     * @enum Precision
     * @brief Enum denoting the precision of operations that are supported to be done on the matrix.
     *
     */
    enum Precision {
        SINGLE = 0,
        DOUBLE = 1,
        MIXED = 2,
    };

    /**
     * @enum FloatPoint
     * @brief Enum denoting the floating point arithmetic of the matrix.
     *
     */
    enum FloatPoint : int {
        EXAGEOSTAT_BYTE = 0,
        EXAGEOSTAT_INTEGER = 1,
        EXAGEOSTAT_REAL_FLOAT = 2,
        EXAGEOSTAT_REAL_DOUBLE = 3,
        EXAGEOSTAT_COMPLEX_FLOAT = 4,
        EXAGEOSTAT_COMPLEX_DOUBLE = 5,
    };

    /**
     * @enum UpperLower
     * @brief Enum denoting the Upper/Lower part
     *
     */
    enum UpperLower : int {
        EXAGEOSTAT_UPPER = 121, /**< Use lower triangle of A */
        EXAGEOSTAT_LOWER = 122, /**< Use upper triangle of A */
        EXAGEOSTAT_UPPER_LOWER = 123  /**< Use the full A */
    };

    /**
     * @enum CopyDirection
     * @brief Enum denoting the copy descriptors flow
     *
     */
    enum CopyDirection : int {
        CHAMELEON_TO_HICMA = 0,
        HICMA_TO_CHAMELEON = 1
    };

    /**
     * @var availableKernels
     * @brief Set denoting the available kernels supported in matrix generation.
     * @details This set is updated automatically to add new kernels.
     * The set is initialized with a lambda function that iterates through a directory
     * and extracts the kernel names from the filenames. It also adds lowercase versions
     * of the kernel names with underscores before each capital letter.
     * @return set of all available kernels names
     *
     */
    const static std::set<std::string> availableKernels = []() {
        // This set stores the kernel names.
        std::set<std::string> kernelNames;
        // This string stores the directory path where the kernel files are located.
        const std::string directoryPath = KERNELS_PATH;
        // This loop iterates through all the files in the directory and extracts the kernel names.
        for (const auto &entry: std::filesystem::directory_iterator(directoryPath)) {
            // This checks if the current entry is a regular file.
            if (entry.is_regular_file()) {
                // This string stores the filename of the current entry.
                const std::string filename = entry.path().filename().string();
                // This string stores the file extension of the current entry.
                const std::string extension = std::filesystem::path(filename).extension().string();
                // This string stores the kernel name extracted from the filename.
                const std::string kernelName = filename.substr(0, filename.size() - extension.size());
                // This adds the kernel name to the kernelNames set.
                kernelNames.insert(kernelName);
                // This blo/ck of code converts the kernel name to lowercase     and adds underscores before each capital letter.
                std::string lowercaseName;
                for (std::size_t i = 0; i < kernelName.size(); ++i) {
                    if (std::isupper(kernelName[i])) {
                        // Avoid adding _ in the beginning of the name.
                        if (i != 0) {
                            lowercaseName += '_';
                        }
                        lowercaseName += static_cast<char>(std::tolower(kernelName[i]));
                    } else {
                        lowercaseName += kernelName[i];
                    }
                }
                // This adds the lowercase kernel name to the kernelNames set.
                kernelNames.insert(lowercaseName);
            }
        }
        return kernelNames;
    }();
}//namespace exageostat

#endif //EXAGEOSTATCPP_DEFINITIONS_HPP
