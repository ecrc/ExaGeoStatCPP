
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file PredictionHelpers.cpp
 * @brief Contains the implementation of the PredictionHelpers class.
 * @version 1.0.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-06-08
**/

#include <algorithm>

#include <prediction/PredictionHelpers.hpp>

using namespace exageostat::prediction;
using namespace exageostat::configurations;
using namespace exageostat::dataunits;

template<typename T>
void PredictionHelpers<T>::PickRandomPoints(Configurations &aConfigurations, ExaGeoStatData<T> &aData, T *apZObs,
                                            T *apZActual, T *apZ, Locations<T> &aMissLocation,
                                            Locations<T> &aObsLocation) {

    int i;
    int j;
    int N = aConfigurations.GetProblemSize();
    int z_miss_number = aConfigurations.GetUnknownObservationsNb();
    int z_obs_number = aConfigurations.CalculateZObsNumber();
    auto l = new Locations<T>(N, aData.GetLocations()->GetDimension());

    int p = aConfigurations.GetP();
    bool is_shuffle = true;

    if (aConfigurations.GetIsMLOEMMOM() && aConfigurations.GetP() == 2) {
        p = 1;
        N /= 2;
        is_shuffle = false;
    }

    // Create an array of pointers using 'new'
    T **Z_parts = new T *[p];
    // Allocate memory for each row of the 2D array
    for (i = 0; i < p; i++) {
        Z_parts[i] = new T[N / p];
    }

    if (p > 1) {
        // Partition Z into p parts
        for (i = 0; i < p; i++) {
            for (j = 0; j < N; j += p) {
                Z_parts[i][j] = apZ[i + j];
            }
        }
    }

    for (i = 0; i < N / p; i++) {
        l->GetLocationX()[i] = aData.GetLocations()->GetLocationX()[i];
        l->GetLocationY()[i] = aData.GetLocations()->GetLocationY()[i];
    }

    if (is_shuffle) {
        if (p == 1) {
            Shuffle(apZ, *l, N);
        } else if (p == 2) {
            Shuffle(Z_parts[0], Z_parts[1], *l, N);
        } else if (p == 3) {
            Shuffle(Z_parts[0], Z_parts[1], Z_parts[2], *l, N);
        }
    }

    j = 0;
    for (i = 0; i < z_miss_number; i++) {
        if (p == 1) {
            apZActual[i] = apZ[i];
        } else if (p == 2) {
            apZActual[j] = Z_parts[0][i];
            apZActual[j + 1] = Z_parts[1][i];
            j += 2;
        } else if (p == 3) {
            apZActual[j] = Z_parts[0][i];
            apZActual[j + 1] = Z_parts[1][i];
            apZActual[j + 2] = Z_parts[2][i];
            j += 3;
        }
    }

    j = 0;
    for (i = 0; i < z_obs_number; i++) {
        if (p == 1) {
            apZObs[i] = apZ[z_miss_number + i];
        } else if (p == 2) {
            apZObs[j] = Z_parts[0][z_miss_number + i];
            apZObs[j + 1] = Z_parts[1][z_miss_number + i];
            j += 2;
        } else if (p == 3) {
            apZObs[j] = Z_parts[0][z_miss_number + i];
            apZObs[j + 1] = Z_parts[1][z_miss_number + i];
            apZObs[j + 2] = Z_parts[2][z_miss_number + i];
            j += 3;
        }
    }

    for (i = 0; i < z_miss_number; i++) {
        aMissLocation.GetLocationX()[i] = l->GetLocationX()[i];
        aMissLocation.GetLocationY()[i] = l->GetLocationY()[i];
    }

    for (i = 0; i < z_obs_number; i++) {
        aObsLocation.GetLocationX()[i] = l->GetLocationX()[z_miss_number + i];
        aObsLocation.GetLocationY()[i] = l->GetLocationY()[z_miss_number + i];
    }

    if (p == 1) {
        SortInplace(z_obs_number, aObsLocation, apZObs);
        SortInplace(z_miss_number, aMissLocation, apZActual);
    }

    for (i = 0; i < p; i++) {
        delete[] Z_parts[i];
    }
    delete[] Z_parts;
    delete l;
}

template<typename T>
void PredictionHelpers<T>::Shuffle(T *apArray, Locations<T> &aLocations, int aSize) {
    if (aSize > 1) {
        size_t i;

        for (i = 0; i < aSize - 1; i++) {
            size_t j = i + rand() / (RAND_MAX / (aSize - i) + 1);

            T t = apArray[j];
            apArray[j] = apArray[i];
            apArray[i] = t;

            T x_temp = aLocations.GetLocationX()[j];
            aLocations.GetLocationX()[j] = aLocations.GetLocationX()[i];
            aLocations.GetLocationX()[i] = x_temp;

            T y_temp = aLocations.GetLocationY()[j];
            aLocations.GetLocationY()[j] = aLocations.GetLocationY()[i];
            aLocations.GetLocationY()[i] = y_temp;

        }
    }
}

template<typename T>
void PredictionHelpers<T>::Shuffle(T *apArray1, T *apArray2, Locations<T> &aLocations, int aSize) {
    if (aSize > 1) {
        size_t i;

        for (i = 0; i < aSize - 1; i++) {
            size_t j = i + rand() / (RAND_MAX / (aSize - i) + 1);

            T t1 = apArray1[j];
            apArray1[j] = apArray1[i];
            apArray1[i] = t1;

            T t2 = apArray2[j];
            apArray2[j] = apArray2[i];
            apArray2[i] = t2;

            T x_temp = aLocations.GetLocationX()[j];
            aLocations.GetLocationX()[j] = aLocations.GetLocationX()[i];
            aLocations.GetLocationX()[i] = x_temp;

            T y_temp = aLocations.GetLocationY()[j];
            aLocations.GetLocationY()[j] = aLocations.GetLocationY()[i];
            aLocations.GetLocationY()[i] = y_temp;

        }
    }

}

template<typename T>
void PredictionHelpers<T>::Shuffle(T *apArray1, T *apArray2, T *apArray3, Locations<T> &aLocations,
                                   int aSize) {
    if (aSize > 1) {
        size_t i;

        for (i = 0; i < aSize - 1; i++) {
            size_t j = i + rand() / (RAND_MAX / (aSize - i) + 1);

            T t1 = apArray1[j];
            apArray1[j] = apArray1[i];
            apArray1[i] = t1;

            T t2 = apArray2[j];
            apArray2[j] = apArray2[i];
            apArray2[i] = t2;

            T t3 = apArray3[j];
            apArray3[j] = apArray3[i];
            apArray3[i] = t3;

            T x_temp = aLocations.GetLocationX()[j];
            aLocations.GetLocationX()[j] = aLocations.GetLocationX()[i];
            aLocations.GetLocationX()[i] = x_temp;

            T y_temp = aLocations.GetLocationY()[j];
            aLocations.GetLocationY()[j] = aLocations.GetLocationY()[i];
            aLocations.GetLocationY()[i] = y_temp;

        }
    }

}

template<typename T>
void PredictionHelpers<T>::SortArray(uint32_t *aData, int aCount) {

    std::sort(aData, aData + aCount);
}

template<typename T>
int PredictionHelpers<T>::SortInplace(int aN, Locations<T> &aLocations, T *apZ) {

    int i;
    int j;//new_j, tmp_j;
    int count = aN;
    int ndim = 2;
    T *point = new T[ndim * count];
    T *ptr1;
    T *minmax; // min is stored in lower part, max is stored in upper part

    for (i = 0; i < count; i++) {
        point[i] = aLocations.GetLocationX()[i];
        point[i + count] = aLocations.GetLocationY()[i];
    }

    minmax = new T[2 * count];

    for (i = 0; i < ndim; i++) {
        ptr1 = point + i * count; // i-th dimension
        minmax[i] = ptr1[0];
        minmax[i + ndim] = minmax[i];
        for (j = 1; j < count; j++) {
            if (minmax[i] > ptr1[j])
                minmax[i] = ptr1[j];
            else if (minmax[i + ndim] < ptr1[j])
                minmax[i + ndim] = ptr1[j];
        }
    }

    // Now minmax[0:ndim] and minmax[ndim:2*ndim] store minimal and maximal
    // values of coordinates
    auto uint_point = new uint32_t[ndim * count];
    uint32_t *uint_ptr1;
    T min, range;

    for (i = 0; i < ndim; i++) {
        uint_ptr1 = uint_point + i * count;
        ptr1 = point + i * count;
        min = minmax[i];
        range = minmax[i + ndim] - min;
        for (j = 0; j < count; j++)
            uint_ptr1[j] = (ptr1[j] - min) / range * UINT32_MAX;
    }

    delete[] minmax;

    // Now uint_ptr1 contains initial coordinates, rescaled to range
    // [0, UINT32_MAX] and converted to uint32_t type to use special radix sort
    // Prepare indexes to store sort order
    auto order = new int[count];
    for (j = 0; j < count; j++)
        order[j] = j;
    SortArray(uint_point, count);

    auto new_point = new T[ndim * count];
    auto new_z = new T[count];

    for (j = 0; j < count; j++) {
        for (i = 0; i < ndim; i++)
            new_point[count * i + j] = point[count * i + order[j]];
        new_z[j] = apZ[order[j]];
    }

    for (i = 0; i < count; i++) {
        aLocations.GetLocationX()[i] = new_point[i];
        aLocations.GetLocationY()[i] = new_point[i + count];
        apZ[i] = new_z[i];
    }

    delete[] new_point;
    delete[] point;
    delete[] new_z;
    delete[] uint_point;
    delete[] order;

    return 0;
}
