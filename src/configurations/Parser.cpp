
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
* @file Parser.cpp
* @version 1.1.0
* @brief Contains the implementation of the Parser class and its member functions for configuration parsing.
* @details Provides static methods to parse command-line arguments and JSON configuration files, as well as utility functions for string transformations.
* @author Mahmoud ElKarargy
* @date 2024-12-11
**/

#include <fstream>
#include <nlohmann/json.hpp>

#include <configurations/Parser.hpp>

using namespace std;
using namespace exageostat::configurations::parser;

void Parser::ParseCLI(const int &aArgC, char **apArgV, unordered_map<string, any> &apConfigurationMap) {

    string example_name = apArgV[0];
    // Remove the './'
    example_name.erase(0, 2);
    string argument;
    string argument_name;
    string argument_value;
    int equal_sign_Idx;

    for (int i = 1; i < aArgC; ++i) {
        argument = apArgV[i];
        argument = argument.substr(2);
        equal_sign_Idx = static_cast<int>(argument.find('='));
        argument_name = argument.substr(0, equal_sign_Idx);

        string converted_name = ProcessKeyString(argument_name);

        if (equal_sign_Idx != string::npos) {
            // Store value as a string
            argument_value = argument.substr(equal_sign_Idx + 1);
            apConfigurationMap[converted_name] = argument_value;
        } else {
            // For flags, store "true" as a string
            apConfigurationMap[converted_name] = string("true");
        }
    }
}

void Parser::ParseJSON(const string &aJSONFile, unordered_map<string, any> &apConfigurationMap) {

    const auto &config_path = aJSONFile;
    ifstream json_file_stream(config_path);
    if (!json_file_stream.is_open()) {
        throw runtime_error("Could not open JSON configuration file.");
    }

    nlohmann::json json_config;
    json_file_stream >> json_config;

    for (auto &[key, value]: json_config.items()) {
        auto converted_name = ProcessKeyString(key);
        apConfigurationMap[converted_name] = value.get<string>();

    }
}

string Parser::ProcessKeyString(const string &arg) {

    // Process the string to remove dashes/underscores and handle capitalization
    std::string result;
    for (char ch : arg) {
        // Ignore hyphens, underscores, and other non-alphanumeric characters
        if (std::isalnum(ch)) { // Keep only alphanumeric characters
            result += static_cast<char>(std::tolower(ch)); // Convert to lowercase
        }
    }
    return result;
}
