
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file DiskWriter.cpp
 * @brief Contains the implementation of the DiskWriter class for writing data to disk.
 * @version 1.0.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-06-08
**/

#include <string>
#include <fstream>
#include <filesystem>

#include <helpers/DiskWriter.hpp>

using namespace std;

using namespace exageostat::helpers;
using namespace exageostat::dataunits;

template<typename T>
void DiskWriter<T>::WriteVectorsToDisk(const T &aMatrixPointer, const int &aProblemSize, const int &aP,
                                       std::string &aLoggerPath, Locations<T> &aLocations) {

    // Determine the path for storing the output files
    if (aLoggerPath.empty()) {
        aLoggerPath = LOG_PATH;
    } else {
        if (aLoggerPath.back() == '/') {
            aLoggerPath += "synthetic_ds";
        } else {
            aLoggerPath += "/synthetic_ds";
        }
    }
    // Create a new directory if it does not already exist
    bool created = false;
    if (!filesystem::exists(aLoggerPath)) {
        try {
            created = filesystem::create_directory(aLoggerPath);
        } catch (const filesystem::filesystem_error &e) {
            throw runtime_error("Error creating directory: " + aLoggerPath);
        }
    } else {
        created = true;
    }

    // Check if the directory was created successfully
    if (!created) {
        throw runtime_error("Error creating directory: " + aLoggerPath);
        return;
    }

    // Determine the names of the output files
    size_t i = 1, j;
    std::ofstream p_file_synthetic, p_file_log;
    std::string n_file_synthetic = aLoggerPath + "/SYN_" + std::to_string(aProblemSize / aP) + "_";
    std::string n_file_log = aLoggerPath + "/log_" + std::to_string(aProblemSize / aP) + "_";
    std::string temp = n_file_log + std::to_string(i);

    // Check if log file exists
    while (std::filesystem::exists(temp)) {
        i++;
        temp = n_file_log + std::to_string(i);
    }

    n_file_synthetic += std::to_string(i);
    p_file_synthetic.open(n_file_synthetic);

    for (j = 0, i = 0; i < aProblemSize / aP; i++) {
        if (aLocations.GetLocationZ() == nullptr) {
            //2 Dimensions
            if (aP == 1) {
                p_file_synthetic << std::setprecision(15) << aLocations.GetLocationX()[i] << ','
                                 << aLocations.GetLocationY()[i] << "," << std::setprecision(15) << (&aMatrixPointer)[i]
                                 << '\n';
            } else if (aP == 2) {
                p_file_synthetic << std::setprecision(15) << aLocations.GetLocationX()[i] << ','
                                 << aLocations.GetLocationY()[i] << "," << std::setprecision(15) << (&aMatrixPointer)[j]
                                 << std::setprecision(15) << (&aMatrixPointer)[j + 1] << '\n';
                j += 2;
            } else if (aP == 3) {
                p_file_synthetic << std::setprecision(15) << aLocations.GetLocationX()[i] << ','
                                 << aLocations.GetLocationY()[i] << "," << std::setprecision(15) << (&aMatrixPointer)[j]
                                 << std::setprecision(15) << (&aMatrixPointer)[j + 1] << std::setprecision(15)
                                 << (&aMatrixPointer)[j + 2] << '\n';
                j += 3;
            }
        } else {
            //3 Dimensions
            if (aP == 1) {
                p_file_synthetic << std::setprecision(15) << aLocations.GetLocationX()[i] << ','
                                 << aLocations.GetLocationY()[i] << ',' << aLocations.GetLocationZ()[i] << ","
                                 << std::setprecision(15) << (&aMatrixPointer)[i] << '\n';
            } else if (aP == 2) {
                p_file_synthetic << std::setprecision(15) << aLocations.GetLocationX()[i] << ','
                                 << aLocations.GetLocationY()[i] << ',' << aLocations.GetLocationZ()[i] << ","
                                 << std::setprecision(15) << (&aMatrixPointer)[j] << std::setprecision(15)
                                 << (&aMatrixPointer)[j + 1] << '\n';
                j += 2;
            } else if (aP == 2) {
                p_file_synthetic << std::setprecision(15) << aLocations.GetLocationX()[i] << ','
                                 << aLocations.GetLocationY()[i] << ',' << aLocations.GetLocationZ()[i] << ","
                                 << std::setprecision(15) << (&aMatrixPointer)[j] << std::setprecision(15)
                                 << (&aMatrixPointer)[j + 1] << std::setprecision(15) << (&aMatrixPointer)[j + 2]
                                 << '\n';
                j += 3;
            }

        }
    }
    p_file_synthetic.close();
}
