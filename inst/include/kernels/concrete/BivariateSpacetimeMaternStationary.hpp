//
//// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
//// Copyright (C) 2023 by Brightskies inc,
//// All rights reserved.
//// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).
//
///**
// * @file BivariateSpacetimeMaternStationary.hpp
// *
// * @version 1.0.0
// * @author Sameh Abdulah
// * @date 2023-04-14
//**/
//
//#ifndef EXAGEOSTATCPP_BIVARIATESPACETIMEMATERNSTATIONARY_HPP
//#define EXAGEOSTATCPP_BIVARIATESPACETIMEMATERNSTATIONARY_HPP
//
//#include <kernels/Kernel.hpp>
//
//namespace exageostat {
//    namespace kernels {
//
//        /**
//         * @class BivariateSpacetimeMaternStationary
//         */
//        class BivariateSpacetimeMaternStationary : public Kernel {
//        public:
//
//            /**
//             * @brief Generates a covariance matrix using a set of locations and kernel parameters.
//             * @param[in] apMatrixA The output covariance matrix.
//             * @param[in] aRowsNumber The number of rows in the output matrix.
//             * @param[in] aColumnsNumber The number of columns in the output matrix.
//             * @param[in] aRowOffset The row offset for the input locations.
//             * @param[in] aColumnOffset The column offset for the input locations.
//             * @param[in] apLocation1 The set of input locations 1.
//             * @param[in] apLocation2 The set of input locations 2.
//             * @param[in] apLocalTheta An array of kernel parameters.
//             * @param [in] aDistanceMetric Distance metric to be used (1 = Euclidean, 2 = Manhattan, 3 = Minkowski).
//             */
//            void GenerateCovarianceMatrix(double *apMatrixA, int aRowsNumber, int aColumnsNumber, int aRowOffset,
//                                          int aColumnOffset, dataunits::Locations *apLocation1,
//                                          dataunits::Locations *apLocation2, dataunits::Locations *apLocation3,
//                                          double *apLocalTheta, int aDistanceMetric) override;
//
////            void InitializeMatrixParameters() override;
//
//        };
//    }//namespace Kernels
//}//namespace exageostat
//
//#endif //EXAGEOSTATCPP_BIVARIATESPACETIMEMATERNSTATIONARY_HPP
