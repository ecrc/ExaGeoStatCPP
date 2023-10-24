
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ChameleonImplementationDST.hpp
 * @brief This file contains the declaration of ChameleonImplementationDST class.
 * @details ChameleonImplementationDST is a concrete implementation of LinearAlgebraMethods class for diagonal super tile matrices.
 * @version 1.0.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-03-26
**/

#ifndef EXAGEOSTATCPP_CHAMELEONIMPLEMENTATIONDST_HPP
#define EXAGEOSTATCPP_CHAMELEONIMPLEMENTATIONDST_HPP

#include <linear-algebra-solvers/concrete/chameleon/ChameleonImplementation.hpp>

namespace exageostat::linearAlgebra::diagonalSuperTile {

    /**
     * @brief ChameleonImplementationDST is a concrete implementation of LinearAlgebraMethods class for diagonal super tile matrices.
     * @tparam T Data Type: float or double
     *
     */
    template<typename T>
    class ChameleonImplementationDST : public ChameleonImplementation<T> {

    public:

        /**
         * @brief Default constructor.
         */
        explicit ChameleonImplementationDST() = default;

        /**
         * @brief Virtual destructor to allow calls to the correct concrete destructor.
         */
        ~ChameleonImplementationDST() override = default;

        /**
         * @brief Computes the Cholesky factorization of a symmetric positive definite or Symmetric positive definite matrix.
         * @copydoc LinearAlgebraMethods::ExaGeoStatPotrfTile()
         */
        void
        ExaGeoStatPotrfTile(const common::UpperLower &aUpperLower, void *apA, int aBand, void *apCD, void *apCrk,
                            const int &aMaxRank, const int &aAcc) override;

        /**
         * @brief Computes the parallel Cholesky factorization of a symmetric positive definite diagonal super tile matrix.
         * @param[in] aUpperLower Whether upper or lower part of the matrix A
         * @param[in] apA Symmetric matrix A
         * @param[in] aBand diagonal thickness.
         * @param[in] apSequence The sequence structure to associate in the options.
         * @param[in] apRequest The request structure to associate in the options.
         * @return successful exit.
         */
        void
        ExaGeoStatParallelPotrfDiag(const common::UpperLower &aUpperLower, void *apA, int aBand, void *apSequence,
                                    void *apRequest);

        /**
         * @brief Computes the Cholesky factorization of a symmetric positive definite diagonal super tile matrix.
         * @param[in] aUpperLower Whether upper or lower part of the matrix A
         * @param[in] apA Symmetric matrix A
         * @param[in] aBand diagonal thickness.
         * @param[in] apSequence The sequence structure to associate in the options.
         * @param[in] apRequest The request structure to associate in the options.
         * @return successful exit.
         */
        int ExaGeoStatPotrfTileAsync(const common::UpperLower &aUpperLower, void *apA, int aBand, void *apSequence,
                                     void *apRequest);
    };

    /**
    * @brief Instantiates the chameleon DST class for float and double types.
    * @tparam T Data Type: float or double
    *
    */
    EXAGEOSTAT_INSTANTIATE_CLASS(ChameleonImplementationDST)

}//namespace exageostat

#endif //EXAGEOSTATCPP_CHAMELEONIMPLEMENTATIONDST_HPP