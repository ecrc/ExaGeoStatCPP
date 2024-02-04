
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file HicmaImplementation.hpp
 * @brief This file contains the declaration of HicmaImplementation class.
 * @details HicmaImplementation is a concrete implementation of LinearAlgebraMethods class for tile low-rank matrices.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-03-26
**/

#ifndef EXAGEOSTATCPP_HICMAIMPLEMENTATION_HPP
#define EXAGEOSTATCPP_HICMAIMPLEMENTATION_HPP

#include <linear-algebra-solvers/LinearAlgebraMethods.hpp>

namespace exageostat::linearAlgebra::tileLowRank {

    /**
     * @brief HicmaImplementation is a concrete implementation of LinearAlgebraMethods class for tile low-rank matrices.
     * @tparam T Data Type: float or double
     *
     */
    template<typename T>
    class HicmaImplementation : public LinearAlgebraMethods<T> {
    public:

        /**
         * @brief Default constructor.
         */
        explicit HicmaImplementation() = default;

        /**
         * @brief Virtual destructor to allow calls to the correct concrete destructor.
         */
        ~HicmaImplementation() override = default;

        /**
         * @brief Set the modeling descriptors for HiCMA implementation.
         * @param[in,out] aData Reference to the ExaGeoStatData object.
         * @param[in] aConfigurations Reference to the Configurations object.
         * @param[in] aP the P value of the kernel multiplied by time slot.
         */
        void SetModelingDescriptors(std::unique_ptr<dataunits::ExaGeoStatData<T>> &aData,
                                    Configurations &aConfigurations, const int &aP);

        /**
         * @brief Calculates the log likelihood value of a given value theta.
         * @copydoc LinearAlgebraMethods::ExaGeoStatMLETile()
        */
        T ExaGeoStatMLETile(const ExaGeoStatHardware &apHardware,
                            std::unique_ptr<dataunits::ExaGeoStatData<T>> &aData,
                            Configurations &aConfigurations, const double *theta,
                            T *apMeasurementsMatrix, const kernels::Kernel<T> &aKernel) override;

        /**
         * @brief Copies a matrix in the tile layout from source to destination
         * @copydoc LinearAlgebraMethods::ExaGeoStatLapackCopyTile()
         */
        void ExaGeoStatLapackCopyTile(const common::UpperLower &aUpperLower, void *apA, void *apB) override;

        /**
         * @brief Initialize the runtime option structure for HiCMA
         * @copydoc LinearAlgebraMethods::ExaGeoStatOptionsInit()
         */
        void
        ExaGeoStatOptionsInit(void *apOptions, void *apContext, void *apSequence, void *apRequest) override;

        /**
         * @brief Submit the release of the workspaces associated to the options structure.
         * @copydoc LinearAlgebraMethods::ExaGeoStatOptionsFree()
         */
        void ExaGeoStatOptionsFree(void *apOptions) override;

        /**
         * @brief Finalize the runtime option structure for HiCMA.
         * @copydoc LinearAlgebraMethods::ExaGeoStatOptionsFinalize()
         */
        void ExaGeoStatOptionsFinalize(void *apOptions, void *apContext) override;

        /**
         * @brief Wait for the completion of a sequence.
         * @copydoc LinearAlgebraMethods::ExaGeoStatSequenceWait()
         */
        void ExaGeoStatSequenceWait(void *apSequence) override;

        /**
         * @brief Create HiCMA Sequence.
         * @copydoc LinearAlgebraMethods::ExaGeoStatCreateSequence()
         */
        void
        ExaGeoStatCreateSequence(void *apSequence) override;

        /**
         * @brief Computes the Cholesky factorization of a symmetric positive definite or Symmetric positive definite matrix.
         * @copydoc LinearAlgebraMethods::ExaGeoStatPotrfTile()
         */
        void
        ExaGeoStatPotrfTile(const common::UpperLower &aUpperLower, void *apA, int aBand, void *apCD, void *apCrk,
                            const int &aMaxRank, const int &aAcc) override;

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
        int ExaGeoStatMeasureDetTileAsync(void *apDescA, void *apSequence, void *apRequest, void *apDescDet) override;

        /**
         * @brief Get the pointer to the data or the runtime handler associated to the piece of data (m, n) in desc.
         * @copydoc LinearAlgebraMethods::ExaGeoStatDataGetAddr()
         */
        void *ExaGeoStatDataGetAddr(void *apA, int aAm, int aAn) override;

        /**
        * @brief Calculate the loglikelihood of non-Gaussian MLE.
        * @copydoc LinearAlgebraMethods::ExaGeoStatNonGaussianLogLikeTileAsync()
        */
        int ExaGeoStatNonGaussianLogLikeTileAsync(void *apDescZ, void *apDescSum, const T *apTheta,
                                                  void *apSequence, void *apRequest) override;

        /**
        * @brief Calculate the loglikelihood of non-Gaussian MLE.
        * @copydoc LinearAlgebraMethods::ExaGeoStatNonGaussianTransformTileAsync()
        */
        int ExaGeoStatNonGaussianTransformTileAsync(void *apDescZ, void *apDescFlag, const T *apTheta,
                                                    void *apSequence, void *apRequest) override;
    };

    /**
    * @brief Instantiates the Hicma TLR class for float and double types.
    * @tparam T Data Type: float or double
    *
    */
    EXAGEOSTAT_INSTANTIATE_CLASS(HicmaImplementation)

}//namespace exageostat

#endif //EXAGEOSTATCPP_HICMAIMPLEMENTATION_HPP