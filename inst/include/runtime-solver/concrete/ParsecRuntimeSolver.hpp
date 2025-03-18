
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ParsecRuntimeSolver.hpp
 * @brief This file contains the declaration of ParsecRuntimeSolver class.
 * @details ParsecRuntimeSolver is a concrete implementation of the RuntimeSolverMethods class.
 * @version 2.0.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @author Qinglei Cao
 * @date 2024-11-04
**/

#ifndef EXAGEOSTATCPP_PARSECRUNTIMESOLVER_HPP
#define EXAGEOSTATCPP_PARSECRUNTIMESOLVER_HPP

#include <runtime-solver/RuntimeSolverMethods.hpp>

namespace exageostat::runtimesolver {

    /**
     * @brief ParsecRuntimeSolver is a concrete implementation of RuntimeSolverMethods class for dense or diagonal-super tile matrices.
     * @tparam T Data Type: float or double
     *
     */
    template<typename T>
    class ParsecRuntimeSolver : public RuntimeSolverMethods<T> {
    public:

        /**
         * @brief Calculates the log likelihood value of a given value theta.
         * @copydoc RuntimeSolverMethods::ModelingOperations()
         *
         */
        T ModelingOperations(std::unique_ptr <ExaGeoStatData<T>> &aData, configurations::Configurations &aConfigurations,
                             T *apMeasurementsMatrix, const kernels::Kernel <T> &aKernel) override;

        /**
         * @brief API function for data modeling.
         * @param aTheta A vector containing the model parameters (theta) to evaluate.
         * @param aGrad A vector to store the computed gradient of the log likelihood with respect to theta.
         * @param apInfo Pointer to additional context information, which may include configuration,
         *               data descriptors, and kernel settings needed for the modeling process.
         * @return The computed log likelihood value.
         */
        static double DataModelingAPI(const std::vector<double> &aTheta, std::vector<double> &aGrad, void *apInfo);

        /**
         * @brief Performs Maximum Likelihood Estimation (MLE) on a tiled matrix using PaRSEC.
         * @param apConfiguration The configuration settings for the data modeling process.
         * @param theta Pointer to an array of doubles representing the current model parameters.
         * @param grad Reference to a vector of doubles to store the gradient computed during modeling.
         * @param data Unique pointer to an ExaGeoStatData instance that holds the matrix data and descriptors.
         * @param aKernel The kernel object containing the statistical model details.
         * @return The log likelihood value computed.
         */
        static double ParsecDmleTile(configurations::Configurations &apConfiguration,const double* theta,
             std::vector<double> & grad, std::unique_ptr<ExaGeoStatData<T>> &data,const kernels::Kernel<T> &aKernel);

        
        /**
         * @brief Performs a SYRK (symmetric rank-k update) operation on the matrix.
         * @param[in,out] aData Descriptor Data object to be populated with descriptors and data.
         *
         */
        static void ExaGeoStatSYRK(std::unique_ptr<ExaGeoStatData<T>> &aData);

        /**
         * @brief Performs TLR Cholesky operation on the matrix.
         * @param[in,out] aData Descriptor Data object to be populated with descriptors and data.
         *
         */
        static void ExaGeoStatTLRCholesky(std::unique_ptr<ExaGeoStatData<T>> &aData);

        /**
         * @brief Calculates norm.
         * @param[in] aConfigurations Configurations object containing relevant settings.
         * @param[in,out] aData Descriptor Data object to be populated with descriptors and data.
         *
         */
        static double ExaGeoStatNorm(configurations::Configurations &aConfigurations, std::unique_ptr<ExaGeoStatData<T>> &aData);

         /**
         * @brief Calculates the Mean Squared Error (MSE).
         * @param[in] aConfigurations Reference to Configurations object containing needed parameters.
         * @param[out] aData Reference to an ExaGeoStatData<T> object that contains matrix to be analyzed.
         * @return the calculated MSE.
         *
         */
        static double CalculateMSE(configurations::Configurations &aConfigurations,
                            std::unique_ptr<ExaGeoStatData<T>> &aData);

    };

    /**
    * @brief Instantiates the Parsec Runtime Solver class for float and double types.
    * @tparam T Data Type: float or double
    *
    */
    EXAGEOSTAT_INSTANTIATE_CLASS(ParsecRuntimeSolver)
}//namespace exageostat

#endif //EXAGEOSTATCPP_PARSECRUNTIMESOLVER_HPP