
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file Validator.cpp
 * @brief This file defines the Validator class which validates configuration string input.
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @date 2024-12-11
**/

#include <nlohmann/json.hpp>

#include <utilities/Logger.hpp>
#include <configurations/Validator.hpp>

using namespace std;
using namespace exageostat::common;
using namespace exageostat::configurations::validator;

const unordered_map<string, function<any(const string &)>> Validator::mCheckersMap = {
        {"Theta",          [](const string &value) -> any { return CheckThetaValue(value); }},
        {"Tolerance",      [](const string &value) -> any { return CheckToleranceValue(value); }},
        {"FileHandler",    [](const string &value) -> any { return CheckLogFileValue(value); }},
        {"Bool",           [](const string &value) -> any { return CheckBoolValue(value); }},
        {"Numerical",      [](const string &value) -> any { return CheckNumericalValue(value); }},
        {"Dimension",      [](const string &value) -> any { return CheckDimensionValue(value); }},
        {"DistanceMetric", [](const string &value) -> any { return CheckDistanceMetricValue(value); }},
        {"Kernel",         [](const string &value) -> any { return CheckKernelValue(value); }},
        {"Verbose",        [](const string &value) -> any { return CheckVerboseValue(value); }},
        {"Precision",      [](const string &value) -> any { return CheckPrecisionValue(value); }},
        {"Computation",    [](const string &value) -> any { return CheckComputationValue(value); }}
};

const std::unordered_map<std::string, std::string> Validator::mArgumentToCategoryMap = {
        {"n",                 "Numerical"},
        {"kernel",            "Kernel"},
        {"p",                 "Numerical"},
        {"q",                 "Numerical"},
        {"timeslot",          "Numerical"},
        {"computation",       "Computation"},
        {"precision",         "Precision"},
        {"cores",             "Numerical"},
        {"gpus",              "Numerical"},
        {"dts",               "Numerical"},
        {"lts",               "Numerical"},
        {"band",              "Numerical"},
        {"maxrank",           "Numerical"},
        {"observationnumber", "Numerical"},
        {"observationsfile",  "Path"},
        {"filelogname",       "Path"},
        {"filelogpath",       "FileHandler"},
        {"seed",              "Numerical"},
        {"logpath",           "Path"},
        {"initialtheta",      "Theta"},
        {"ooc",               "Bool"},
        {"approximationmode", "Numerical"},
        {"mspe",              "Bool"},
        {"fisher",            "Bool"},
        {"idw",               "Bool"},
        {"log",               "Bool"},
        {"lb",                "Theta"},
        {"ub",                "Theta"},
        {"estimatedtheta",    "Theta"},
        {"startingtheta",     "Theta"},
        {"isnongaussian",     "Bool"},
        {"verbose",           "Verbose"},
        {"banddensedp",       "Numerical"},
        {"banddense",         "Numerical"},
        {"banddensehp",       "Numerical"},
        {"objectsnumber",     "Numerical"},
        {"adaptivedecision",  "Numerical"},
        {"adddiagonal",       "Numerical"},
        {"filetimeslot",      "Numerical"},
        {"filenumber",        "Numerical"},
        {"hnb",               "Numerical"},
        {"genmaxrank",        "Numerical"},
        {"compmaxrank",       "Numerical"},
        {"autoband",          "Numerical"},
        {"banddensesp",       "Numerical"},
        {"bandlowrankdp",     "Numerical"},
        {"enableinverse",     "Bool"},
        {"mpiio",             "Bool"},
        {"mloemmom",          "Bool"},
        {"dimension",         "Dimension"},
        {"stagezero",         "Bool"},
        {"issynthetic",       "Bool"},
        {"datapath",          "Path"},
        {"resultspath",       "Path"},
        {"lat",               "Numerical"},
        {"lon",               "Numerical"},
        {"netcdfdatapath",        "Path"},
        {"forcingdatapath",        "Path"},
        {"startyear",         "Numerical"},
        {"endyear",           "Numerical"},
        {"isclimateemulator", "Bool"},
        {"recoveryfile",      "Path"},
        {"distancemetric",    "DistanceMetric"},
        {"maxmleiterations",  "Numerical"},
        {"accuracy",          "Numerical"},
        {"tolerance",         "Tolerance"},
        {"zmiss",             "Numerical"},
        {"missingpattern",    "Theta"},
};

Computation Validator::CheckComputationValue(const string &aValue) {

    if (aValue != "exact" and aValue != "Exact" and aValue != "Dense" and aValue != "dense" and
        aValue != "diag_approx" and aValue != "diagonal_approx" and aValue != "lr_approx" and aValue != "tlr" and
        aValue != "TLR") {
        throw range_error("Invalid value for Computation. Please use Exact, diagonal_approx or TLR.");
    }
    if (aValue == "exact" or aValue == "Exact" or aValue == "Dense" or aValue == "dense") {
        return EXACT_DENSE;
    } else if (aValue == "diag_approx" or aValue == "diagonal_approx") {
        return DIAGONAL_APPROX;
    }
    return TILE_LOW_RANK;
}

Precision Validator::CheckPrecisionValue(const string &aValue) {

    if (aValue != "single" and aValue != "Single" and aValue != "double" and aValue != "Double" and aValue != "mix" and
        aValue != "Mix" and aValue != "Mixed" and aValue != "mixed") {
        throw range_error("Invalid value for Computation. Please use Single, Double or Mixed.");
    }
    if (aValue == "single" or aValue == "Single") {
        return SINGLE;
    } else if (aValue == "double" or aValue == "Double") {
        return DOUBLE;
    }
    return MIXED;
}

string Validator::CheckKernelValue(const string &aKernel) {

    // Check if the kernel name exists in the availableKernels set.
    if (availableKernels.count(aKernel) <= 0) {
        throw range_error("Invalid value for Kernel. Please check manual.");
    }
    // Check if the string is already in CamelCase format
    if (IsCamelCase(aKernel)) {
        return aKernel;
    }
    string str = aKernel;
    // Replace underscores with spaces and split the string into words
    replace(str.begin(), str.end(), '_', ' ');
    istringstream iss(str);
    string word, result;
    while (iss >> word) {
        // Capitalize the first letter of each word and append it to the result
        word[0] = static_cast<char>(toupper(word[0]));
        result += word;
    }
    return result;
}

bool Validator::IsCamelCase(const string &aString) {

    // If the string contains an underscore, it is not in CamelCase format
    if (aString.find('_') != string::npos) {
        return false;
    }
    // If the string starts with a lowercase letter, it is not in CamelCase format
    if (islower(aString[0])) {
        return false;
    }
    // If none of the above conditions hold, the string is in CamelCase format
    return true;
}

Verbose Validator::CheckVerboseValue(const string &aVerbosity) {

    if (aVerbosity == "quiet" || aVerbosity == "Quiet") {
        Configurations::SetVerbosity(Verbose::QUIET_MODE);
        return Verbose::QUIET_MODE;
    } else if (aVerbosity == "standard" || aVerbosity == "Standard") {
        Configurations::SetVerbosity(Verbose::STANDARD_MODE);
        return Verbose::STANDARD_MODE;
    } else if (aVerbosity == "detailed" || aVerbosity == "Detailed" || aVerbosity == "detail") {
        Configurations::SetVerbosity(Verbose::DETAILED_MODE);
        return Verbose::DETAILED_MODE;
    } else {
        LOGGER("Error: " << aVerbosity << " is not valid ")
        throw range_error("Invalid value. Please use verbose or standard values only.");
    }
}


int Validator::CheckNumericalValue(const string &aValue) {

    int numerical_value;
    try {
        numerical_value = stoi(aValue);
    }
    catch (...) {
        throw range_error("Invalid value. Please use Numerical values only.");
    }

    if (numerical_value < 0) {
        throw range_error("Invalid value. Please use positive values");
    }
    return numerical_value;
}


vector<double> Validator::CheckThetaValue(const string &aInputValues) {

    // Count the number of values in the string
    int num_values = 1;
    for (char input_value_char: aInputValues) {
        if (input_value_char == ':') {
            num_values++;
        }
    }
    // Allocate memory for the array of doubles
    vector<double> theta;

    if (aInputValues.empty()) {
        return theta;
    }

    // Split the string into tokens using strtok()
    const char *delim = ":";
    char *token = strtok((char *) aInputValues.c_str(), delim);
    int i = 0;
    while (token != nullptr) {
        // Check if the token is a valid double or "?"
        if (!strcmp(token, "?")) {
            theta.push_back(-1);
        } else {
            try {
                theta.push_back(stod(token));
            }
            catch (...) {
                LOGGER("Error: " << token << " is not a valid double or '?' ")
                throw range_error("Invalid value. Please use Numerical values only.");
            }
        }

        // Get the next token
        token = strtok(nullptr, delim);
        i++;
    }

    // Check if the number of values in the array is correct
    if (i != num_values) {
        throw range_error(
                "Error: the number of values in the input string is invalid, please use this example format as a reference 1:?:0.1");
    }

    return theta;
}

Dimension Validator::CheckDimensionValue(const string &aDimension) {

    if (aDimension != "2D" and aDimension != "2d" and aDimension != "3D" and aDimension != "3d" and
        aDimension != "st" and aDimension != "ST") {
        throw range_error("Invalid value for Dimension. Please use 2D, 3D or ST.");
    }
    if (aDimension == "2D" or aDimension == "2d") {
        return Dimension2D;
    } else if (aDimension == "3D" or aDimension == "3d") {
        return Dimension3D;
    }
    return DimensionST;
}

bool Validator::CheckBoolValue(const string &aBooleanValue) {

    // Convert input to lowercase for case-insensitive comparison
    std::string lower_case_value = aBooleanValue;
    std::transform(lower_case_value.begin(), lower_case_value.end(), lower_case_value.begin(), ::tolower);

    // Check if the value matches "true" in any case
    if (lower_case_value == "true")
        return true;
    else
        return false;
}


DistanceMetric Validator::CheckDistanceMetricValue(const string &aDistanceMetric) {

    if (aDistanceMetric == "eg" || aDistanceMetric == "EG" || aDistanceMetric == "euclidean") {
        return EUCLIDEAN_DISTANCE;
    } else if (aDistanceMetric == "gcd" || aDistanceMetric == "GCD" || aDistanceMetric == "great_circle") {
        return GREAT_CIRCLE_DISTANCE;
    } else {
        throw range_error("Invalid value. Please use eg or gcd values only.");
    }
}

double Validator::CheckToleranceValue(const string &aTolerancePower) {

    auto val = Validator::CheckNumericalValue(aTolerancePower);
    return pow(10, -1 * val);
}

FILE *Validator::CheckLogFileValue(const std::string &apFileLogPath) {

    if (apFileLogPath.empty()) {
        throw std::invalid_argument("Path cannot be empty.");
    }

    auto file_handler = fopen(apFileLogPath.c_str(), "w");

    if (file_handler == nullptr) {
        throw std::domain_error("Failed to open file at path " + apFileLogPath + ".");
    }
    return file_handler;
}

void Validator::Validate(unordered_map<string, any> &apArgsMap) {

    for (auto &[key, value]: apArgsMap) {
        // Find the appropriate checker based on key
        auto category_map_key = mArgumentToCategoryMap.find(key);
        if (category_map_key == mArgumentToCategoryMap.end()) { // throw exception if no category
            throw invalid_argument("No category for the following argument name: " + key);
        } else if (category_map_key->second != "Path") {
            auto it = mCheckersMap.find(category_map_key->second);
            if (it != mCheckersMap.end()) {
                if (value.type() == typeid(std::string)) {
                    // Apply the checker function
                    auto checker_function = it->second;
                    auto checked_value = checker_function(any_cast<string>(value));
                    apArgsMap[key] = checked_value;
                }
            } else {
                throw invalid_argument("Unknown key in the args map: " + key);
            }
        }
    }
}
