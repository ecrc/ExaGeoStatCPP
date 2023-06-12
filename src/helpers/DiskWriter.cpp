
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
#include <helpers/DiskWriter.hpp>
#include <sys/stat.h>
#include <fstream>
#include <cstring>
#include <vector>
#include <unistd.h>

using namespace exageostat::helpers;
using namespace std;
using namespace exageostat::dataunits;
//template<typename T>
//void DiskWriter<T>::CreateDirectory(const std::string& path) {
//    struct stat st = { 0 };
//    if (stat(path.c_str(), &st) == -1)
//        mkdir(path.c_str(), 0700);
//}
//
//
//template<typename T>
//int DiskWriter<T>::DoesFileExist(const char* filename) {
//    ifstream infile(filename);
//    return infile.good();
//}

template<typename T>
void DiskWriter<T>::WriteVectorsToDisk(T *apMatrixPointer, const int *apProblemSize, const int *apP, std::string *apLoggerPath, Locations *apLocations) {

    int i = 1;
    FILE* pFileZ = nullptr;
    FILE* pFileZ2 = nullptr;
    FILE* pFileZ3 = nullptr;
    FILE* pFileXY = nullptr;
    struct stat st = {0};

    char nFileZ[100];
    char nFileZ2[100];
    char nFileZ3[100];
    char temp[100];
    char nFileXY[100];
    char nFileLog[100];
    int p = 1;

    // Create new directory if not exist
    if (stat("./synthetic_ds", &st) == -1) {
        mkdir("./synthetic_ds", 0700);
    }
    
    snprintf(nFileZ, 100, "%s%d%s%s%s", "./synthetic_ds/Z1_", *apProblemSize / *apP, "_");
    snprintf(nFileZ2, 100, "%s%d%s%s%s", "./synthetic_ds/Z2_", *apProblemSize / *apP, "_");
    snprintf(nFileZ3, 100, "%s%d%s%s%s", "./synthetic_ds/Z3_", *apProblemSize / *apP, "_");
    snprintf(nFileXY, 100, "%s%d%s%s%s", "./synthetic_ds/LOC_", *apProblemSize / *apP, "_");
    snprintf(nFileLog, 100, "%s%d%s%s%s", "./synthetic_ds/log_", *apProblemSize / *apP, "_");

    snprintf(temp, 100, "%s%d", nFileLog, i);

    // Check if log file exists
    while (access(temp, F_OK) != -1) {
        i++;
        snprintf(temp, 100, "%s%d", nFileLog, i);
    }

    sprintf(temp, "%d", i);
    strcat(nFileZ, temp);
    strcat(nFileXY, temp);
    strcat(nFileZ2, temp);
    strcat(nFileZ3, temp);
    strcat(nFileLog, temp);

    pFileZ = fopen(nFileZ, "w+");
    pFileXY = fopen(nFileXY, "w+");

    if (p == 1) {
        for (i = 0; i < *apProblemSize; i++) {
            fprintf(pFileZ, "%f\n", apMatrixPointer[i]);
        }
    } else if (p == 2) {
        pFileZ2 = fopen(nFileZ2, "w+");
        for (i = 0; i < *apProblemSize; i += 2) {
            fprintf(pFileZ, "%f\n", apMatrixPointer[i]);
            fprintf(pFileZ2, "%f\n", apMatrixPointer[i + 1]);
        }
        fclose(pFileZ2);
    } else if (p == 3) {
        pFileZ2 = fopen(nFileZ2, "w+");
        pFileZ3 = fopen(nFileZ3, "w+");
        for (i = 0; i < *apProblemSize; i += 3) {
            fprintf(pFileZ, "%f\n", apMatrixPointer[i]);
            fprintf(pFileZ2, "%f\n", apMatrixPointer[i + 1]);
            fprintf(pFileZ3, "%f\n", apMatrixPointer[i+ 2]);
        }
        fclose(pFileZ2);
        fclose(pFileZ3);
    }

    for (i = 0; i < *apProblemSize / p; i++) {
        if (apLocations->GetLocationZ() == nullptr) {
            fprintf(pFileXY, "%f,%f\n", apLocations->GetLocationX(), apLocations->GetLocationY());
        } else {
            fprintf(pFileXY, "%f,%f,%f\n", apLocations->GetLocationX(), apLocations->GetLocationY(), apLocations->GetLocationZ());
        }
    }

    fclose(pFileZ);
    fclose(pFileXY);
}
