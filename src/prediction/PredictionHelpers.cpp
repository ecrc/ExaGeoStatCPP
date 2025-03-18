
// Copyright (c) 2017-2024 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file PredictionHelpers.cpp
 * @brief Contains the implementation of the PredictionHelpers class.
 * @version 1.1.0
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
void PredictionHelpers<T>::PickRandomPoints(Configurations &aConfigurations,
                                            std::unique_ptr<ExaGeoStatData<T>> &aData, T *apZObs,
                                            T *apZActual, T *apZ, Locations<T> &aMissLocation,
                                            Locations<T> &aObsLocation, const int &aP) {

    int i;
    int j;
    int full_problem_size = aConfigurations.GetProblemSize() * aP;
    int z_miss_number = aConfigurations.GetUnknownObservationsNb();
    int z_obs_number = aConfigurations.CalculateZObsNumber();
    bool is_shuffle = true;
    int p = aP;
    auto l = new Locations<T>(aConfigurations.GetProblemSize(), aData->GetLocations()->GetDimension());

    if (aConfigurations.GetIsMLOEMMOM() && aP == 2) {
        p = 1;
        full_problem_size /= 2;
        is_shuffle = false;
    }

    // Create an array of pointers using 'new'
    T **Z_parts = new T *[p];
    // Allocate memory for each row of the 2D array
    for (i = 0; i < p; i++) {
        Z_parts[i] = new T[full_problem_size / p];
    }

    if (p > 1) {
        // Partition Z into p parts
        for (i = 0; i < p; i++) {
            int m = 0;
            for (j = 0; j < full_problem_size; j += p) {
                Z_parts[i][m] = apZ[i + j];
                m++;
            }
        }
    }

    for (i = 0; i < full_problem_size / p; i++) {
        l->GetLocationX()[i] = aData->GetLocations()->GetLocationX()[i];
        l->GetLocationY()[i] = aData->GetLocations()->GetLocationY()[i];
        if (aConfigurations.GetDimension() != common::Dimension2D) {
            l->GetLocationZ()[i] = aData->GetLocations()->GetLocationZ()[i];
        }
    }

    if (is_shuffle) {
        if (p == 1) {
            Shuffle(apZ, *l, full_problem_size);
        } else if (p == 2) {
            Shuffle(Z_parts[0], Z_parts[1], *l, full_problem_size / p);
        } else if (p == 3) {
            Shuffle(Z_parts[0], Z_parts[1], Z_parts[2], *l, full_problem_size / p);
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
        if (aConfigurations.GetDimension() != common::Dimension2D) {
            aMissLocation.GetLocationZ()[i] = l->GetLocationZ()[i];
        }
    }

    for (i = 0; i < z_obs_number; i++) {
        aObsLocation.GetLocationX()[i] = l->GetLocationX()[z_miss_number + i];
        aObsLocation.GetLocationY()[i] = l->GetLocationY()[z_miss_number + i];
        if (aConfigurations.GetDimension() != common::Dimension2D) {
            aObsLocation.GetLocationZ()[i] = l->GetLocationZ()[z_miss_number + i];
        }
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

            if (aLocations.GetDimension() != common::Dimension2D) {
                T z_temp = aLocations.GetLocationZ()[j];
                aLocations.GetLocationZ()[j] = aLocations.GetLocationZ()[i];
                aLocations.GetLocationZ()[i] = z_temp;
            }

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

            if (aLocations.GetDimension() != common::Dimension2D) {
                T z_temp = aLocations.GetLocationZ()[j];
                aLocations.GetLocationZ()[j] = aLocations.GetLocationZ()[i];
                aLocations.GetLocationZ()[i] = z_temp;
            }

        }
    }

}

template<typename T>
void PredictionHelpers<T>::Shuffle(T *apArray1, T *apArray2, T *apArray3, Locations<T> &aLocations, int aSize) {
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

            if (aLocations.GetDimension() != common::Dimension2D) {
                T z_temp = aLocations.GetLocationZ()[j];
                aLocations.GetLocationZ()[j] = aLocations.GetLocationZ()[i];
                aLocations.GetLocationZ()[i] = z_temp;
            }
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
    int dimension_number = 2;
    T *point = new T[dimension_number * count];
    T *ptr1;
    T *minmax; // min is stored in lower part, max is stored in upper part

    for (i = 0; i < count; i++) {
        point[i] = aLocations.GetLocationX()[i];
        point[i + count] = aLocations.GetLocationY()[i];
    }

    minmax = new T[2 * count];

    for (i = 0; i < dimension_number; i++) {
        ptr1 = point + i * count; // i-th dimension
        minmax[i] = ptr1[0];
        minmax[i + dimension_number] = minmax[i];
        for (j = 1; j < count; j++) {
            if (minmax[i] > ptr1[j])
                minmax[i] = ptr1[j];
            else if (minmax[i + dimension_number] < ptr1[j])
                minmax[i + dimension_number] = ptr1[j];
        }
    }

    // Now minmax[0:dimension_number] and minmax[dimension_number:2*dimension_number] store minimal and maximal
    // values of coordinates
    auto uint_point = new uint32_t[dimension_number * count];
    uint32_t *uint_ptr1;
    T min, range;

    for (i = 0; i < dimension_number; i++) {
        uint_ptr1 = uint_point + i * count;
        ptr1 = point + i * count;
        min = minmax[i];
        range = minmax[i + dimension_number] - min;
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

    auto new_point = new T[dimension_number * count];
    auto new_z = new T[count];

    for (j = 0; j < count; j++) {
        for (i = 0; i < dimension_number; i++)
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

#if !DEFAULT_RUNTIME
template<typename T>
void PredictionHelpers<T>::RadixSortRecursive(uint32_t *apData, int aCount, int aNDim,
    int *apOrder, int *apTmpOrder, int aSDim, int aSBit, int aLo, int aHi) {
    int aI, aLoLast = aLo, aHiLast = aHi;
    uint32_t *apSData = apData + aSDim * aCount;
    uint32_t aCheck = 1 << aSBit;
    for(aI = aLo; aI <= aHi; aI++)
    {
        if((apSData[apOrder[aI]] & aCheck) == 0)
        {
            apTmpOrder[aLoLast] = apOrder[aI];
            aLoLast++;
        }
        else
        {
            apTmpOrder[aHiLast] = apOrder[aI];
            aHiLast--;
        }
    }
    for(aI = aLo; aI <= aHi; aI++)
        apOrder[aI] = apTmpOrder[aI];
    if(aSDim > 0)
    {
        if(aLoLast - aLo > 1)
        PredictionHelpers<T>::RadixSortRecursive(apData, aCount, aNDim, apOrder, apTmpOrder, aSDim - 1,
                aSBit, aLo, aLoLast - 1);
        if(aHi - aHiLast > 1)
        PredictionHelpers<T>::RadixSortRecursive(apData, aCount, aNDim, apOrder, apTmpOrder, aSDim - 1,
                aSBit, aHiLast + 1, aHi);
    }
    else if(aSBit > 0)
    {
        if(aLoLast - aLo > 1)
        PredictionHelpers<T>::RadixSortRecursive(apData, aCount, aNDim, apOrder, apTmpOrder, aNDim - 1,
                aSBit - 1, aLo, aLoLast - 1);
        if(aHi - aHiLast > 1)
        PredictionHelpers<T>::RadixSortRecursive(apData, aCount, aNDim, apOrder, apTmpOrder, aNDim - 1,
                aSBit - 1, aHiLast + 1, aHi);
    }
}

template<typename T>
int PredictionHelpers<T>::RadixSort(uint32_t *apData, int aCount, int aNDim, int *apOrder)
{
    int *apTmpOrder = (int *) malloc(aCount * sizeof(int));
    PredictionHelpers<T>::RadixSortRecursive(apData, aCount, aNDim, apOrder, apTmpOrder, aNDim - 1, 31, 0, aCount - 1);
    free(apTmpOrder);
    return 0;
}

template<typename T>
int PredictionHelpers<T>::LocationsZSortInPlace(int aN, location *apLocations, double *apZ)
{
    int aI, aJ;
    int aCount = aN;
    int aNDim = 2;
    int aInfo;
    double *apPoint = (double *) malloc(aNDim * aCount * sizeof(double));
    double *apPtr1;
    double *apMinMax; // min values stored in first half, max values in second half

    for(aI = 0; aI < aCount; aI++)
    {
        apPoint[aI]        = apLocations->x[aI];
        apPoint[aI + aCount] = apLocations->y[aI];
    }

    apMinMax = (double *) malloc(2 * aCount * sizeof(double));
    for(aI = 0; aI < aNDim; aI++)
    {
        apPtr1 = apPoint + aI * aCount;
        apMinMax[aI] = apPtr1[0];
        apMinMax[aI + aNDim] = apMinMax[aI];
        for(aJ = 1; aJ < aCount; aJ++)
        {
            if(apMinMax[aI] > apPtr1[aJ])
                apMinMax[aI] = apPtr1[aJ];
            else if (apMinMax[aI + aNDim] < apPtr1[aJ])
                apMinMax[aI + aNDim] = apPtr1[aJ];
        }
    }
    // apMinMax[0:aNDim] stores minima, apMinMax[aNDim:2*aNDim] stores maxima.
    uint32_t *apUintPoint = (uint32_t *) malloc(aNDim * aCount * sizeof(uint32_t));
    uint32_t *apUintPtr1;
    double aMin, aRange;
    for(aI = 0; aI < aNDim; aI++)
    {
        apUintPtr1 = apUintPoint + aI * aCount;
        apPtr1 = apPoint + aI * aCount;
        aMin = apMinMax[aI];
        aRange = apMinMax[aI + aNDim] - aMin;
        for(aJ = 0; aJ < aCount; aJ++)
            apUintPtr1[aJ] = (apPtr1[aJ] - aMin) / aRange * UINT32_MAX;
    }
    free(apMinMax);
    // Prepare indices for sorting.
    int *apOrder = (int *) malloc(aCount * sizeof(int));
    for(aJ = 0; aJ < aCount; aJ++)
        apOrder[aJ] = aJ;
    aInfo = PredictionHelpers<T>::RadixSort(apUintPoint, aCount, aNDim, apOrder);
    if(aInfo != 0)
    {
        free(apUintPoint);
        free(apOrder);
        return aInfo;
    }
    double *apNewPoint = (double *) malloc(aNDim * aCount * sizeof(double));
    double *apNewZ     = (double *) malloc(aCount * sizeof(double));
    for(aJ = 0; aJ < aCount; aJ++)
    {
        for(aI = 0; aI < aNDim; aI++)
            apNewPoint[aCount * aI + aJ] = apPoint[aCount * aI + apOrder[aJ]];
        apNewZ[aJ] = apZ[apOrder[aJ]];
    }
    for(aI = 0; aI < aCount; aI++)
    {
        apLocations->x[aI] = apNewPoint[aI];
        apLocations->y[aI] = apNewPoint[aI + aCount];
        apZ[aI] = apNewZ[aI];
    }
    free(apNewPoint);
    free(apPoint);
    free(apNewZ);
    free(apUintPoint);
    free(apOrder);
    return 0;
}
#endif