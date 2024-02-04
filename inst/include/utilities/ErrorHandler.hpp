// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ErrorHandler.hpp
 * @version 1.1.0
 * @brief Provides error handling functionalities.
 * @details Defines macros and functions for handling errors and warnings.
 * @author Mahmoud ElKarargy
 * @author David Helmy
 * @date 2024-01-20
**/

#ifndef EXAGEOSTATCPP_ERRORHANDLER_HPP
#define EXAGEOSTATCPP_ERRORHANDLER_HPP

#include <utilities/Printer.hpp>

#ifdef USE_CUDA
#include <cuda_runtime.h>
#endif

/**
 * @brief EXAGEOSTAT API Exceptions Macro to use for Errors and Warnings.
 */
 #define API_EXCEPTION(MESSAGE, ERROR_TYPE) \
    APIException(MESSAGE, ERROR_TYPE)

#ifdef USE_CUDA
/**
 * @brief Useful macro wrapper for all cuda API calls to ensure correct returns,
 * and error throwing on failures.
 */
#define GPU_ERROR_CHECK(ans) { APIException::AssertGPU((ans), __FILE__, __LINE__); }
#endif

/**
 * @brief Enumeration for error types.
 */
enum ErrorType : int {
    RUNTIME_ERROR = 0,
    RANGE_ERROR = 1,
    INVALID_ARGUMENT_ERROR = 2,
    WARNING = 3,
};

/**
 * @class APIException
 * @brief Custom exception class for handling API errors and warnings.
 */
class APIException : public std::exception{

public:

    /**
     * @brief Constructor for APIException.
     * @param[in] aMessage The error or warning message.
     * @param[in] aErrorCode The error type.
     */
    APIException(const std::string &aMessage, const ErrorType &aErrorCode) {

        if (aErrorCode != WARNING) {
            APIException::ThrowError(aMessage, aErrorCode);
        } else {
            APIException::ThrowWarning(aMessage);
        }
    }

    /**
     * @brief Destructor for APIException.
     */
    ~APIException() override = default;

#ifdef USE_CUDA
    /**
     * @brief Function to assert the return code of a CUDA API call and ensure it completed successfully.
     * @param[in] aCode The code returned from the CUDA API call.
     * @param[in] aFile The name of the file that the assertion was called from.
     * @param[in] aLine The line number in the file that the assertion was called from.
     */
    inline static void AssertGPU(cudaError_t aCode, const char *aFile, int aLine)
    {
        if (aCode != cudaSuccess)
        {
#ifdef USING_R
            std::string s="GPU Assert:  "+std::string(cudaGetErrorString(aCode));
            Rcpp::stop(s);
#else
            char s[200];
            sprintf((char*)s,"GPU Assert: %s %s %d\n", cudaGetErrorString(aCode), aFile, aLine);
            throw std::invalid_argument(s);
#endif
        }
    }
#endif


private:

    /**
     * @brief Helper function to throw error based on error type.
     * @param[in] aString The error message.
     * @param[in] aErrorCode The error type.
     */
    static void ThrowError(const std::string &aString, const ErrorType &aErrorCode) {
#ifdef USING_R
//        Rcpp::stop(aString);
#else
        if (aErrorCode == RUNTIME_ERROR) {
            throw std::runtime_error(aString);
        } else if (aErrorCode == INVALID_ARGUMENT_ERROR) {
            throw std::invalid_argument(aString);
        } else if (aErrorCode == RANGE_ERROR) {
            throw std::range_error(aString);
        }
#endif
    }

    /**
     * @brief Helper function to throw warning.
     * @param[in] aString The warning message.
     */
    static void ThrowWarning(const std::string &aString) {
#ifdef USING_R
//        Rcpp::warning(aString);
#endif
    }
};

#endif //EXAGEOSTATCPP_ERRORHANDLER_HPP