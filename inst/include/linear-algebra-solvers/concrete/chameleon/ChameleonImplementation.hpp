
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ChameleonImplementation.hpp
 * @brief This file contains the declaration of ChameleonImplementation class.
 * @details ChameleonImplementation is a concrete implementation of the LinearAlgebraMethods class for the common functionality implementation shared between dense and diagonal-super tile matrices.
 * @version 1.0.0
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
         */
        T ExaGeoStatMLETile(const hardware::ExaGeoStatHardware &aHardware, dataunits::ExaGeoStatData<T> &aData,
                            configurations::Configurations &aConfigurations, const double *theta,
                            T *apMeasurementsMatrix, const kernels::Kernel<T> &aKernel) override;

        /**
         * @brief Copy Lapack matrix to Descriptor Matrix
         * @copydoc LinearAlgebraMethods::ExaGeoStatLap2Desc()
         */
        void ExaGeoStatLap2Desc(T *apA, const int &aLDA, void *apDescA, const common::UpperLower &aUpperLower) override;

        /**
         * @brief Copies a matrix in the tile layout from source to destination
         * @copydoc LinearAlgebraMethods::ExaGeoStatLapackCopyTile()
         */
        void ExaGeoStatLapackCopyTile(const common::UpperLower &aUpperLower, void *apA, void *apB) override;

        /**
         * @brief Get the pointer to the data or the runtime handler associated to the piece of data (m, n) in desc.
         * @copydoc LinearAlgebraMethods::ExaGeoStatDataGetAddr()
         */
        void *ExaGeoStatDataGetAddr(void *apA, int aAm, int aAn) override;

        /**
         * @brief  Solves one of the matrix equations op( A )*X = alpha*B, or X*op( A ) = alpha*B.
         * @copydoc LinearAlgebraMethods::ExaGeoStatTrsmTile()
         */
        void ExaGeoStatTrsmTile(const common::Side &aSide, const common::UpperLower &aUpperLower,
                                const common::Trans &aTrans, const common::Diag &aDiag, const T &aAlpha, void *apA,
                                void *apCD, void *apCrk, void *apZ, const int &aMaxRank) override;

        /**
         * @brief Calculate determinant for triangular matrix.
         * @copydoc LinearAlgebraMethods::ExaGeoStatMeasureDetTileAsync()
         */
        int
        ExaGeoStatMeasureDetTileAsync(void *apDescA, void *apSequence, void *apRequest, void *apDescDet) override;

        /**
         * @brief Wait for the completion of a sequence.
         * @copydoc LinearAlgebraMethods::ExaGeoStatSequenceWait()
         */
        void
        ExaGeoStatSequenceWait(void *apSequence) override;

        /**
         * @brief Create CHAMELEON Sequence.
         * @copydoc LinearAlgebraMethods::ExaGeoStatCreateSequence()
         */
        void
        ExaGeoStatCreateSequence(void *apSequence) override;

        /**
         * @brief Initialize the runtime option structure for CHAMELEON
         * @copydoc LinearAlgebraMethods::ExaGeoStatOptionsInit()
         */
        void
        ExaGeoStatOptionsInit(void *apOptions, void *apContext, void *apSequence, void *apRequest) override;

        /**
         * @brief Submit the release of the workspaces associated to the options structure.
         * @copydoc LinearAlgebraMethods::ExaGeoStatOptionsFree()
         */
        void
        ExaGeoStatOptionsFree(void *apOptions) override;

        /**
         * @brief Finalize the runtime option structure for CHAMELEON.
         * @copydoc LinearAlgebraMethods::ExaGeoStatOptionsFinalize()
         */
        void
        ExaGeoStatOptionsFinalize(void *apOptions, void *apContext) override;

        /**
         * @brief copy Chameleon descriptor to vector float*.
         * @param[in] apDescA Exageostat descriptor A.
         * @param[in] apDescB Exageostat descriptor B.
         * @param[in] apDescC Exageostat descriptor C.
         * @param[in] apSequence Identifies the sequence of function calls that this call belongs to.
         * @param[in] apRequest Identifies this function call (for exception handling purposes).
         * @return Returns 0 for success, error code otherwise.
         *
         */
        int
        ExaGeoStaStrideVectorTileAsync(void *apDescA, void *apDescB, void *apDescC, void *apSequence, void *apRequest);

        /**
         * @brief Copy Chameleon descriptor to vector float*.
         * @param[in] apDescA Exageostat descriptor A.
         * @param[in] apDescB Exageostat descriptor B.
         * @param[in] apDescC Exageostat descriptor C.
         * @param[in] apDescD Exageostat descriptor D.
         * @param[in] apSequence Identifies the sequence of function calls that this call belongs to.
         * @param[in] apRequest Identifies this function call (for exception handling purposes).
         * @return Returns 0 for success, error code otherwise.
         */
        int
        ExaGeoStaStrideVectorTileAsync(void *apDescA, void *apDescB, void *apDescC, void *apDescD, void *apSequence,
                                       void *apRequest);

    };

    EXAGEOSTAT_INSTANTIATE_CLASS(ChameleonImplementation)
}//namespace exageostat

#endif //EXAGEOSTATCPP_CHAMELEONIMPLEMENTATION_HPP