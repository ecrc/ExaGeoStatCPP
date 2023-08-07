
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file DiskWriter.cpp
 * @brief 
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-06-08
**/

#include <cstring>
#include <string>
#include <fstream>

#include <helpers/DiskWriter.hpp>

using namespace exageostat::helpers;
using namespace exageostat::dataunits;

using namespace std;

template<typename T>
void DiskWriter<T>::WriteVectorsToDisk(T *apMatrixPointer, const int *apProblemSize, const int *apP,
                                       std::string aLoggerPath, Locations<T> *apLocations) {

    // Determine the path for storing the output files
    if (aLoggerPath.empty()) {
        aLoggerPath = "./synthetic_ds";
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
            cerr << "Error creating directory: " << e.what() << endl;
        }
    } else {
        created = true;
    }

    // Check if the directory was created successfully
    if (!created) {
        cerr << "Error creating directory: " << aLoggerPath << endl;
        return;
    }

    // Determine the names of the output files
    size_t i = 1;
    std::ofstream p_file_z, p_file_z2, p_file_z3, p_file_xy, p_file_log;
    std::string n_file_z = aLoggerPath + "/Z1_" + std::to_string(*apProblemSize / *apP) + "_";
    std::string n_file_z2 = aLoggerPath + "/Z2_" + std::to_string(*apProblemSize / *apP) + "_";
    std::string n_file_z3 = aLoggerPath + "/Z3_" + std::to_string(*apProblemSize / *apP) + "_";
    std::string n_file_xy = aLoggerPath + "/LOC_" + std::to_string(*apProblemSize / *apP) + "_";
    std::string n_file_log = aLoggerPath + "/log_" + std::to_string(*apProblemSize / *apP) + "_";
    std::string temp = n_file_log + std::to_string(i);

    // Check if log file exists
    while (std::filesystem::exists(temp)) {
        i++;
        temp = n_file_log + std::to_string(i);
    }

    n_file_z += std::to_string(i);
    n_file_xy += std::to_string(i);
    n_file_z2 += std::to_string(i);
    n_file_z3 += std::to_string(i);

    p_file_z.open(n_file_z);
    p_file_xy.open(n_file_xy);

    if (*apP == 1) {
        for (i = 0; i < *apProblemSize; i++) {
            p_file_z << std::setprecision(15) << apMatrixPointer[i] << '\n';
        }
    } else if (*apP == 2) {
        p_file_z2.open(n_file_z2);
        for (i = 0; i < *apProblemSize; i += 2) {
            p_file_z << std::setprecision(15) << apMatrixPointer[i] << '\n';
            p_file_z2 << std::setprecision(15) << apMatrixPointer[i + 1] << '\n';
        }
        p_file_z2.close();
    } else if (*apP == 3) {
        p_file_z2.open(n_file_z2);
        p_file_z3.open(n_file_z3);
        for (i = 0; i < *apProblemSize; i += 3) {
            p_file_z << std::setprecision(15) << apMatrixPointer[i] << '\n';
            p_file_z2 << std::setprecision(15) << apMatrixPointer[i + 1] << '\n';
            p_file_z3 << std::setprecision(15) << apMatrixPointer[i + 2] << '\n';
        }
        p_file_z2.close();
        p_file_z3.close();
    }

    for (i = 0; i < *apProblemSize / *apP; i++) {
        if (apLocations->GetLocationZ() == nullptr) {
            p_file_xy << std::setprecision(15) << apLocations->GetLocationX()[i] << ','
                      << apLocations->GetLocationY()[i] << '\n';
        } else {
            p_file_xy << std::setprecision(15) << apLocations->GetLocationX()[i] << ','
                      << apLocations->GetLocationY()[i] << ',' << apLocations->GetLocationZ()[i] << '\n';
        }
    }

    // Close the output files
    p_file_z.close();
    p_file_z2.close();
    p_file_z3.close();
    p_file_xy.close();
}
