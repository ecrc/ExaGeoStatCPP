
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file SyntheticGenerator.hpp
 * @brief A class for generating synthetic data.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-02-14
**/

#ifndef EXAGEOSTAT_CPP_SYNTHETICGENERATOR_HPP
#define EXAGEOSTAT_CPP_SYNTHETICGENERATOR_HPP

#include <data-generators/DataGenerator.hpp>

namespace exageostat::generators::synthetic {

    /**
     * @class SyntheticGenerator
     * @brief A class for generating synthetic data.
     * @tparam T Data Type: float or double
     * @details This class generates synthetic data for use in testing machine learning models.
     *
     */
    template<typename T>
    class SyntheticGenerator : public DataGenerator<T> {

    public:

        /**
         * @brief Get a pointer to the singleton instance of the SyntheticGenerator class.
         * @return A pointer to the instance of the SyntheticGenerator class.
         *
         */
        static SyntheticGenerator<T> *GetInstance();

        /**
         * @brief Creates the data by synthetically generating it.
         * @copydoc DataGenerator::CreateData()
         *
         */
        std::unique_ptr<ExaGeoStatData<T>>
        CreateData(configurations::Configurations &aConfigurations,
                   exageostat::kernels::Kernel<T> &aKernel) override;

        /**
         * @brief Release the singleton instance of the SyntheticGenerator class.
         * @return void
         *
         */
        static void ReleaseInstance();

    private:
        /**
         * @brief Constructor for the SyntheticGenerator class.
         * @return void
         *
         */
        SyntheticGenerator() = default;

        /**
         * @brief Default destructor.
         *
         */
        ~SyntheticGenerator() override = default;

        /**
         * @brief Pointer to the singleton instance of the SyntheticGenerator class.
         *
         */
        static SyntheticGenerator<T> *mpInstance;

    };

    /**
     * @brief Instantiates the Synthetic Data Generator class for float and double types.
     * @tparam T Data Type: float or double
     *
     */
    EXAGEOSTAT_INSTANTIATE_CLASS(SyntheticGenerator)
} // namespace exageostat

#endif //EXAGEOSTAT_CPP_SYNTHETICGENERATOR_HPP