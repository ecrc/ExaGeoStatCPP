
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file CSVDataGenerator.cpp
 * @brief Implementation of the CSVDataGenerator class
 * @version 1.1.0
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-02-14
**/

#include <fstream>

#include <data-loader/concrete/CSVLoader.hpp>

using namespace std;

using namespace exageostat::dataunits;
using namespace exageostat::kernels;
using namespace exageostat::common;
using namespace exageostat::linearAlgebra;
using namespace exageostat::dataLoader::csv;

template<typename T>
CSVLoader<T> *CSVLoader<T>::GetInstance() {

    if (mpInstance == nullptr) {
        mpInstance = new CSVLoader<T>();
    }
    return mpInstance;
}

template<typename T>
void CSVLoader<T>::ReadData(Configurations &aConfigurations, vector<T> &aMeasurementsMatrix, vector<T> &aXLocations,
                           vector<T> &aYLocations, vector<T> &aZLocations, const int &aP) {

    //Check if the user entered a valid path for the CSV file.
    if (aConfigurations.GetDataPath().empty()) {
        throw runtime_error("Please enter the path to data file.");
    }

    string data_path = aConfigurations.GetDataPath();
    Dimension dimension = aConfigurations.GetDimension();
    ifstream file;

    file.open(data_path, ios::in);

    //Throw error if unable to open file.
    if (!file.is_open()) {
        throw runtime_error("Cannot read locations file: " + data_path);
    }

    //Keep track of number of lines.
    int index = 0;
    string line;
    //Loop through the lines
    while (getline(file, line)) {
        istringstream iss(line);
        string token;

        // Split each line into tokens using ',' as the delimiter
        if (getline(iss, token, ',')) {
            aXLocations.push_back(stod(token));
        }
        if (getline(iss, token, ',')) {
            aYLocations.push_back(stod(token));
        }
        if (getline(iss, token, ',')) {
            //If its a 2D locations' data, the last values of the lines should be the measurement values.
            //If its 3D locations' data, the third value of the line should be the Z coordinate.
            //If its ST location data, the third value of the line should be the Time coordinate
            if (dimension == Dimension2D) {
                aMeasurementsMatrix.push_back(stod(token));
                //if p == 2, the third and fourth values of each line are saved in the Measurements Matrix.
                if (aP == 2) {
                    if (getline(iss, token, ',')) {
                        aMeasurementsMatrix.push_back(stod(token));
                    } else {
                        throw runtime_error(
                                "The data P in the provided file isn't consistent with the Kernel's P.");
                    }
                }
                if (aP == 3) {
                    //if p == 3, the third, fourth, and fifth values of each line are saved in the Measurements Matrix.
                    if (getline(iss, token, ',')) {
                        aMeasurementsMatrix.push_back(stod(token));
                    }
                    if (getline(iss, token, ',')) {
                        aMeasurementsMatrix.push_back(stod(token));
                    } else {
                        throw runtime_error(
                                "The data P in the provided file isn't consistent with the Kernel's P.");
                    }
                }
            } else {
                aZLocations.push_back(stod(token));
            }
        }

        if (getline(iss, token, ',')) {
            if (dimension != Dimension2D) {
                aMeasurementsMatrix.push_back(stod(token));
                //if p == 2, the fourth and fifth values of each line are saved in the Measurements Matrix.
                if (aP == 2) {
                    if (getline(iss, token, ',')) {
                        aMeasurementsMatrix.push_back(stod(token));
                    } else {
                        throw runtime_error(
                                "The data P in the provided file isn't consistent with the Kernel's P.");
                    }
                }
                if (aP == 3) {
                    //if p == 3, the fourth, fifth, and sixth values of each line are saved in the Measurements Matrix.
                    if (getline(iss, token, ',')) {
                        aMeasurementsMatrix.push_back(stod(token));
                    }
                    if (getline(iss, token, ',')) {
                        aMeasurementsMatrix.push_back(stod(token));
                    } else {
                        throw runtime_error(
                                "The data P in the provided file isn't consistent with the Kernel's P.");
                    }
                }
            } else {
                throw runtime_error(
                        "The data dimensions in the provided file isn't consistent with the dimensions input.");
            }
        } else if (dimension != Dimension2D) {
            throw runtime_error("The data dimensions in the provided file isn't consistent with the dimensions input.");
        }
        index++;
    }

    //problem size is equal to the number of the CSV lines * P / aConfigurations.GetTimeSlot().
    aConfigurations.SetProblemSize(index * aP / aConfigurations.GetTimeSlot());

    file.close();
    LOGGER("Data is read from " << data_path << " successfully.")
}

template<typename T>
void CSVLoader<T>::ReleaseInstance() {
    if (mpInstance != nullptr) {
        mpInstance = nullptr;
    }
}

template<typename T> CSVLoader<T> *CSVLoader<T>::mpInstance = nullptr;
