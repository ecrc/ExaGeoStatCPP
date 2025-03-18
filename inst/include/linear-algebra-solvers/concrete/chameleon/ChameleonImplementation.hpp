
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ChameleonImplementation.hpp
 * @brief This file contains the declaration of ChameleonImplementation class.
 * @details ChameleonImplementation is a concrete implementation of the LinearAlgebraMethods class for the common functionality implementation shared between dense and diagonal-super tile matrices.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-03-20
**/

#ifndef EXAGEOSTATCPP_CHAMELEONIMPLEMENTATION_HPP
#define EXAGEOSTATCPP_CHAMELEONIMPLEMENTATION_HPP

#include <linear-algebra-solvers/LinearAlgebraMethods.hpp>

namespace exageostat::linearAlgebra {

    /**
     * @brief ChameleonImplementation is a concrete implementation of LinearAlgebraMethods class for dense or diagonal-super tile matrices.
     * @tparam T Data Type: float or double
     *
     */
    template<typename T>
    class ChameleonImplementation : public LinearAlgebraMethods<T> {
    public:

        /**
         * @brief Calculates the log likelihood value of a given value theta.
         * @copydoc LinearAlgebraMethods::ExaGeoStatMLETile()
         *
         */
        T ExaGeoStatMLETile(std::unique_ptr<ExaGeoStatData<T>> &aData,
                            configurations::Configurations &aConfigurations, const double *theta,
                            T *apMeasurementsMatrix, const kernels::Kernel<T> &aKernel) override;

        /**
         * @brief Copies a matrix in the tile layout from source to destination
         * @copydoc LinearAlgebraMethods::ExaGeoStatLapackCopyTile()
         *
         */
        void ExaGeoStatLapackCopyTile(const common::UpperLower &aUpperLower, void *apA, void *apB) override;

        /**
         * @brief  Solves one of the matrix equations op( A )*X = alpha*B, or X*op( A ) = alpha*B.
         * @copydoc LinearAlgebraMethods::ExaGeoStatTrsmTile()
         *
         */
        void ExaGeoStatTrsmTile(const common::Side &aSide, const common::UpperLower &aUpperLower,
                                const common::Trans &aTrans, const common::Diag &aDiag, const T &aAlpha, void *apA,
                                void *apCD, void *apCrk, void *apZ, const int &aMaxRank) override;

        /**
         * @brief Wait for the completion of a sequence.
         * @copydoc LinearAlgebraMethods::ExaGeoStatSequenceWait()
         *
         */
        void
        ExaGeoStatSequenceWait(void *apSequence) override;

        /**
         * @brief Create CHAMELEON Sequence.
         * @copydoc LinearAlgebraMethods::ExaGeoStatCreateSequence()
         *
         */
        void
        ExaGeoStatCreateSequence(void *apSequence) override;
    };

    /**
    * @brief Instantiates the Chameleon Implementation class for float and double types.
    * @tparam T Data Type: float or double
    *
    */
    EXAGEOSTAT_INSTANTIATE_CLASS(ChameleonImplementation)
}//namespace exageostat

#endif //EXAGEOSTATCPP_CHAMELEONIMPLEMENTATION_HPP