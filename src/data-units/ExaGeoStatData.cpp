/**
 * @file ExaGeoStatData.cpp
 * @brief 
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-07-21
**/
#include <data-units/ExaGeoStatData.hpp>
#include <iostream>

using namespace exageostat::dataunits;
using namespace exageostat::common;

template<typename T>
ExaGeoStatData<T>::ExaGeoStatData(int aSize, Dimension aDimension) {
    auto *locations = new Locations<T>(aSize, aDimension);
    this->mpLocations = locations;
    auto *descriptorData = new DescriptorData<T>;
    this->mpDescriptorData = descriptorData;
}

template<typename T>
ExaGeoStatData<T>::~ExaGeoStatData() {
    delete this->mpLocations;
    delete this->mpDescriptorData;
}

template<typename T>
Locations<T> *ExaGeoStatData<T>::GetLocations() {
    return this->mpLocations;
}

template<typename T>
void ExaGeoStatData<T>::SetLocations(Locations<T> *apLocation) {
    this->mpLocations = apLocation;
    this->mpLocations->SetLocationX(apLocation->GetLocationX());
    this->mpLocations->SetLocationY(apLocation->GetLocationY());
}

template<typename T>
DescriptorData<T> *ExaGeoStatData<T>::GetDescriptorData() {
    return this->mpDescriptorData;
}

template<typename T>
Locations<T> *ExaGeoStatData<T>::CalculateMedianLocations(std::string &aKernelName) {

    auto median_Locations = new Locations<T>(this->mpLocations->GetSize(), this->mpLocations->GetDimension());

    if (aKernelName == "UnivariateMaternNonStationary") {

        T x_min = this->mpLocations->GetLocationX()[0], x_max = this->mpLocations->GetLocationX()[0],
                y_min = this->mpLocations->GetLocationY()[0], y_max = this->mpLocations->GetLocationY()[0],
                z_min, z_max;


        if (this->mpLocations->GetDimension() != common::Dimension2D) {
            z_min = this->mpLocations->GetLocationZ()[0];
            z_max = this->mpLocations->GetLocationZ()[0];
        }

        for (int i = 0; i < this->mpLocations->GetSize(); i++) {
            T x = this->mpLocations->GetLocationX()[i];
            T y = this->mpLocations->GetLocationY()[i];

            x_min = (x < x_min) ? x : x_min;
            x_max = (x > x_max) ? x : x_max;
            y_min = (y < y_min) ? y : y_min;
            y_max = (y > y_max) ? y : y_max;

            if (this->mpLocations->GetDimension() != common::Dimension2D) {
                T z = this->mpLocations->GetLocationX()[i];
                z_min = (z < z_min) ? z : z_min;
                z_max = (z > z_max) ? z : z_max;
            }
        }

        median_Locations->GetLocationX()[0] = x_min + (x_max - x_min) / 2;
        median_Locations->GetLocationY()[0] = y_min + (y_max - y_min) / 2;
        if (this->mpLocations->GetDimension() != common::Dimension2D) {
            median_Locations->GetLocationZ()[0] = z_min + (z_max - z_min) / 2;
        }
    } else {
        median_Locations->GetLocationX()[0] = 0.5;
        median_Locations->GetLocationY()[0] = 0.5;
        if (this->mpLocations->GetDimension() != common::Dimension2D) {
            median_Locations->GetLocationY()[0] = 0.5;
        }
    }
    return median_Locations;
}