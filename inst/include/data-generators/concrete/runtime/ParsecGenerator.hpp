
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file ParsecGenerator.hpp
 * @brief A class for generating synthetic data for Parsec, it implements synthetic generator interface.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2025-02-17
**/

#ifndef EXAGEOSTAT_CPP_PARSECGENERATOR_HPP
#define EXAGEOSTAT_CPP_PARSECGENERATOR_HPP

#include <data-generators/concrete/SyntheticGenerator.hpp>

namespace exageostat::generators::synthetic::parsec {

    /**
     * @class ParsecGenerator
     * @brief A class for generating synthetic data for Parsec.
     * @tparam T Data Type: float or double
     * @details This class generates synthetic data for Parsec runtime to be used in testing machine learning models.
     *
     */
    template<typename T>
    class ParsecGenerator : public SyntheticGenerator<T> {

    public:

        /**
         * @brief Get a pointer to the singleton instance of the ParsecGenerator class.
         * @return A pointer to the instance of the ParsecGenerator class.
         *
         */
        static ParsecGenerator<T> *GetInstance();

        /**
         * @brief Creates the data synthetically.
         * @copydoc DataGenerator::CreateData()
         *
         */
        std::unique_ptr<ExaGeoStatData<T>> 
        CreateSyntheticData(configurations::Configurations &aConfigurations, exageostat::kernels::Kernel<T> &aKernel) override;

        /**
         * @brief Release the singleton instance of the ParsecGenerator class.
         * @return void
         *
         */
        static void ReleaseInstance();

    private:
        /**
         * @brief Constructor for the ParsecGenerator class.
         * @return void
         *
         */
        ParsecGenerator() = default;

        /**
         * @brief Default destructor.
         *
         */
        ~ParsecGenerator() override = default;

        /**
         * @brief Pointer to the singleton instance of the ParsecGenerator class.
         *
         */
        static ParsecGenerator<T> *mpInstance;

    };

    /**
     * @brief Instantiates the Parsec Synthetic Data Generator class for float and double types.
     * @tparam T Data Type: float or double
     *
     */
    EXAGEOSTAT_INSTANTIATE_CLASS(ParsecGenerator)
} // namespace exageostat

#endif //EXAGEOSTAT_CPP_PARSECGENERATOR_HPP