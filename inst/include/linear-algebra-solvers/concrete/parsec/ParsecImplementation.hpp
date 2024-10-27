
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ParsecImplementation.hpp
 * @brief This file contains the declaration of ParsecImplementation class.
 * @details ParsecImplementation is a concrete implementation of the LinearAlgebraMethods class.
 * @version 2.0.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2024-10-15
**/

#ifndef EXAGEOSTATCPP_PARSECIMPLEMENTATION_HPP
#define EXAGEOSTATCPP_PARSECIMPLEMENTATION_HPP

#include <linear-algebra-solvers/LinearAlgebraMethods.hpp>

namespace exageostat::linearAlgebra {

    /**
     * @brief ParsecImplementation is a concrete implementation of LinearAlgebraMethods class.
     * @tparam T Data Type: float or double
     *
     */
    template<typename T>
    class ParsecImplementation : public LinearAlgebraMethods<T> {

        /**
         * @brief Performs a SYRK (symmetric rank-k update) operation on the matrix.
         * @copydoc LinearAlgebraMethods::ExaGeoStatSYRK()
         *
         */
        void ExaGeoStatSYRK(configurations::Configurations &aConfigurations, std::unique_ptr<ExaGeoStatData<T>> &aData) override;

        /**
         * @brief Performs TLR Cholesky operation on the matrix.
         * @copydoc LinearAlgebraMethods::ExaGeoStatTLRCholesky()
         *
         */
        void ExaGeoStatTLRCholesky(configurations::Configurations &aConfigurations, std::unique_ptr<ExaGeoStatData<T>> &aData) override;

        /**
         * @brief Calculates norm.
         * @copydoc LinearAlgebraMethods::ExaGeoStatNorm()
         *
         */
        void ExaGeoStatNorm(configurations::Configurations &aConfigurations,
                            std::unique_ptr<ExaGeoStatData<T>> &aData) override;

        /**
         * @brief Calculates the Mean Squared Error (MSE).
         * @copydoc LinearAlgebraMethods::CalculateMSE()
         *
         */
        double CalculateMSE(configurations::Configurations &aConfigurations,
                            std::unique_ptr<ExaGeoStatData<T>> &aData) override;

        /**
         * @brief The Gateway for the Modeling Operation
         * @copydoc LinearAlgebraMethods::ModelingOperations()
         *
         */
        virtual T ModelingOperations(std::unique_ptr<ExaGeoStatData<T>> &aData,
                                     configurations::Configurations &aConfigurations, const double *apTheta,
                                     T *apMeasurementsMatrix, const kernels::Kernel<T> &aKernel) = 0;

        /**
         * @brief Copies a matrix in the tile layout from source to destination
         * @copydoc LinearAlgebraMethods::ExaGeoStatLapackCopyTile()
         *
         */
        virtual void ExaGeoStatLapackCopyTile(const common::UpperLower &aUpperLower, void *apA, void *apB) = 0;

        /**
         * @brief Wait for the completion of a sequence.
         * @copydoc LinearAlgebraMethods::ExaGeoStatSequenceWait()
         *
         */
        virtual void ExaGeoStatSequenceWait(void *apSequence) = 0;

        /**
         * @brief Create Sequence.
         * @copydoc LinearAlgebraMethods::ExaGeoStatCreateSequence()
         *
         */
        virtual void
        ExaGeoStatCreateSequence(void *apSequence) = 0;

        /**
        * @brief Computes the Cholesky factorization of a symmetric positive definite or Symmetric positive definite matrix.
         * @copydoc LinearAlgebraMethods::ExaGeoStatPotrfTile()
        *
        */
        virtual void
        ExaGeoStatPotrfTile(const common::UpperLower &aUpperLower, void *apA, int aBand, void *apCD, void *apCrk,
                            const int &aMaxRank, const int &aAcc) = 0;

        /**
         * @brief Solves one of the matrix equations op( A )*X = alpha*B, or X*op( A ) = alpha*B.
         * @copydoc LinearAlgebraMethods::ExaGeoStatTrsmTile()
         *
         */
        virtual void ExaGeoStatTrsmTile(const common::Side &aSide, const common::UpperLower &aUpperLower,
                                        const common::Trans &aTrans, const common::Diag &aDiag, const T &aAlpha,
                                        void *apA, void *apCD, void *apCrk, void *apZ, const int &aMaxRank) = 0;

        /**
         * @brief Copy Descriptor Matrix to another Descriptor matrix.
         * @copydoc LinearAlgebraMethods::CopyDescriptors()
         *
         */
        virtual void CopyDescriptors(void *apSourceDesc, void *apDestinationDesc, const int &aSize,
                                     const common::CopyDirection &aDirection) = 0;

    };
    EXAGEOSTAT_INSTANTIATE_CLASS(ParsecImplementation)
}//namespace exageostat

#endif //EXAGEOSTATCPP_PARSECIMPLEMENTATION_HPP