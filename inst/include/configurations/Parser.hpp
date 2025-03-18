
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
* @file Parser.hpp
* @version 1.1.0
* @brief Contains the declaration of the Parser class and its member functions for configuration parsing.
* @details Provides static methods to parse command-line arguments and JSON configuration files, as well as utility functions for string transformations.
* @author Mahmoud ElKarargy
* @date 2024-12-11
**/

#ifndef EXAGEOSTAT_CPP_PARSER_HPP
#define EXAGEOSTAT_CPP_PARSER_HPP

#include <unordered_map>
#include <any>

#include <common/Definitions.hpp>

namespace exageostat::configurations::parser {

    /**
     * @class Parser
     * @brief A class containing static methods for parsing configurations and utility functions.
     */
    class Parser {
    public:

        /**
         * @brief Parses command-line arguments and extracts them into a key-value map.
         * @param[in] aArgC The number of command-line arguments.
         * @param[in] apArgV The array of command-line arguments.
         * @param[in] apConfigurationMap The configuration map to fill.
         * @return void.
         *
         */
        static void ParseCLI(const int &aArgC, char **apArgV, std::unordered_map<std::string, std::any> &apConfigurationMap);

        /**
         * @brief Parses a JSON configuration file and extracts its contents into a key-value map.
         * @param[in] aJSONFile The path to the JSON file.
         * @param[in] apConfigurationMap The configuration map to fill.
         * @return void.
         *
         */
        static void ParseJSON(const std::string &aJSONFile, std::unordered_map<std::string, std::any> &apConfigurationMap);

        /**
         * @brief Converts a given string to all-small characters format.
         * @param[in] arg The input string to convert.
         * @return The converted string.
         */
        static std::string ProcessKeyString(const std::string &arg);

    };
}// namespace exageostat::configurations::parser

#endif // EXAGEOSTAT_CPP_PARSER_HPP
