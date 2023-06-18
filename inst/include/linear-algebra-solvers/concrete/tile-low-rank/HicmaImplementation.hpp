
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file HicmaImplementation.hpp
 * @brief This file contains the declaration of HicmaImplementation class.
 * @details HicmaImplementation is a concrete implementation of LinearAlgebraMethods class for tile low-rank matrices.
 * @version 1.0.0
 * @date 2023-03-26
**/

#ifndef EXAGEOSTATCPP_HICMAIMPLEMENTATION_HPP
#define EXAGEOSTATCPP_HICMAIMPLEMENTATION_HPP

#include <linear-algebra-solvers/LinearAlgebraMethods.hpp>

namespace exageostat {
    namespace linearAlgebra {
        namespace tileLowRank {

            /**
             * @brief HicmaImplementation is a concrete implementation of LinearAlgebraMethods class for tile low-rank matrices.
             * @tparam T Data Type: float or double
             * 
             */
            template<typename T>
            class HicmaImplementation : public LinearAlgebraMethods<T>{
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
                 * @brief Initializes the descriptors necessary for the linear algebra solver.
                 * @copydoc LinearAlgebraMethods::InitiateDescriptors()
                 * 
                 */
                void InitiateDescriptors() override;

                /**
                 * @brief Destroys the descriptors used by the linear algebra solver.
                 * @copydoc LinearAlgebraMethods::DestoryDescriptors()
                 * 
                 */
                void DestoryDescriptors() override;

                /**
                 * @brief Computes the covariance matrix.
                 * @copydoc LinearAlgebraMethods::CovarianceMatrixCodelet()
                 * 
                 */
                void
                CovarianceMatrixCodelet(void *apDescriptor, int &aTriangularPart, dataunits::Locations *apLocation1,
                                        dataunits::Locations *apLocation2,
                                        dataunits::Locations *apLocation3, double *aLocalTheta, int aDistanceMetric,
                                        exageostat::kernels::Kernel *apKernel) override;

                /**
                 * @brief Generates the observations vector.
                 * @copydoc LinearAlgebraMethods::GenerateObservationsVector()
                 * 
                 */
                void GenerateObservationsVector(void *apDescriptor, dataunits::Locations *apLocation1,
                                                dataunits::Locations *apLocation2,
                                                dataunits::Locations *apLocation3, std::vector<double> aLocalTheta,
                                                int aDistanceMetric, exageostat::kernels::Kernel *apKernel) override;

                /**
                 * @brief Initializes the context needed for the Chameleon solver.
                 * @copydoc LinearAlgebraMethods::ExaGeoStatInitContext()
                 * 
                 */
                void ExaGeoStatInitContext(const int &apCoresNumber, const int &apGPUs) override;

                /**
                 * @brief Finalizes the context needed for the Chameleon solver.
                 * @copydoc LinearAlgebraMethods::ExaGeoStatFinalizeContext()
                 * 
                 */
                void ExaGeoStatFinalizeContext() override;

                /**
                 * @brief Copies the descriptor data to a double vector.
                 * @copydoc LinearAlgebraMethods::CopyDescriptorZ()
                 *
                 */
                void CopyDescriptorZ(void *apapDescriptor, double *apDoubleVector) override;


                /**
                 * @brief allocates dense matrix tile.
                 * @copydoc LinearAlgebraMethods::ExageostatAllocateMatrixTile()
                 * 
                 */
                void ExageostatAllocateMatrixTile(void **apDescriptor, bool aIsOOC, T *apMemSpace, int aType2, int aMB,
                                                  int aNB, int aMBxNB, int aLda, int aN, int aSMB, int aSNB, int aM,
                                                  int aN2, int aP, int aQ) override;


            private:
                //// Used context
                static void *apContext;
            };

            /**
            * @brief Instantiates the Hicma TLR class for float and double types.
            * @tparam T Data Type: float or double
            *
            */
            EXAGEOSTAT_INSTANTIATE_CLASS(HicmaImplementation);

        }//namespace tileLowRank
    }//namespace linearAlgebra
}//namespace exageostat

#endif //EXAGEOSTATCPP_HICMAIMPLEMENTATION_HPP