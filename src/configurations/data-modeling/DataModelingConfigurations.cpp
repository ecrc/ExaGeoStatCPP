
/*
 * Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
 * All rights reserved.
 * ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).
 */

/**
 * @file DataModelingConfigurations.cpp
 * @brief 
 * @version 1.0.0
 * @author Sameh Abdulah
 * @date 2023-06-21
**/

#include <configurations/data-modeling/DataModelingConfigurations.hpp>

using namespace std;

using namespace exageostat::configurations::data_modeling;


void DataModelingConfigurations::InitModuleArguments(int aArgC, char **apArgV) {
    InitializeArguments(aArgC, apArgV);
    this->mpMatrixDeterminant = new double;

    this->mDotProduct.resize(3);
    this->mVariance.resize(3);

    this->mDotProduct[0] = new double;
    this->mDotProduct[1] = new double;
    this->mDotProduct[2] = new double;

    *this->mDotProduct[0] = 0;
    *this->mDotProduct[1] = 0;
    *this->mDotProduct[2] = 0;

    this->mVariance[0] = new double;
    this->mVariance[1] = new double;
    this->mVariance[2] = new double;

    cout << "Mproduct Done" <<endl;
    string argument;
    string argument_name;
    string argument_value;
    int equal_sign_Idx;
    // Loop through the arguments that are specific for data generation.
    for (int i = 1; i < aArgC; ++i) {
        argument = apArgV[i];
        equal_sign_Idx = static_cast<int>(argument.find('='));
        argument_name = argument.substr(0, equal_sign_Idx);

        // Check if argument has an equal sign.
        if (equal_sign_Idx != string::npos) {
            argument_value = argument.substr(equal_sign_Idx + 1);

            // Check the argument name and set the corresponding value
            if (argument_name == "--iterations" || argument_name == "--Iterations") {
                SetIterationsValue(CheckNumericalValue(argument_value));
            } else if (argument_name == "--distance_metric" || argument_name == "--distanceMetric"){
                ParseDistanceMetric(argument_value);
            } else if (argument_name == "--data_log" || argument_name == "--dataLog"){
                ParseDataLog(argument_value);
            }
        } else {
//            if (argument_name == "--emptyForNow" || argument_name == "--EmptyForNow") {
//            }
        }
    }
}

DataModelingConfigurations::DataModelingConfigurations(const DataModelingConfigurations &aDataModelingConfigurations) {

    // Call to default copy constructor
    *this = aDataModelingConfigurations;
    InitModuleArguments(this->mArgC, this->mpArgV);
}

DataModelingConfigurations::DataModelingConfigurations(int aArgC, char **apArgV) {

    InitializeArguments(aArgC, apArgV);
    InitModuleArguments(aArgC, apArgV);
}

void DataModelingConfigurations::SetRecoveryFile(string &aRecoveryFile) {
    this->mRecoveryFile = aRecoveryFile;
}

string DataModelingConfigurations::GetRecoveryFile() {
    return this->mRecoveryFile;
}

void DataModelingConfigurations::SetIterationsValue(int aIterationValue) {
    this->mIterationsValue = aIterationValue;
}

int DataModelingConfigurations::GetIterationsValue() const {
    return this->mIterationsValue;
}

double *DataModelingConfigurations::GetDeterminantValue() {
    if (this->mpMatrixDeterminant == nullptr) {
        throw std::runtime_error("Matrix determinant is null");
    }
    return this->mpMatrixDeterminant;
}

std::vector<double *> &DataModelingConfigurations::GetDotProduct() {
    if( this->mDotProduct[0] == nullptr){
        cout<< "get product [0] is null" <<endl;
    }
    return this->mDotProduct;
}

void DataModelingConfigurations::ParseDistanceMetric(const std::string &aDistanceMetric) {
    if (aDistanceMetric == "eg" || aDistanceMetric == "EG"){
        SetDistanceMetric("eg");
    } else if (aDistanceMetric == "gcd" || aDistanceMetric == "GCD"){
        SetDistanceMetric("gcd");
    }else {
        throw range_error("Invalid value. Please use eg or gcd values only.");
    }
}

void DataModelingConfigurations::SetDistanceMetric (const std::string &aDistanceMetric){
        this->mDistanceMetric = aDistanceMetric;
}

std::string
DataModelingConfigurations::GetDistanceMetric (){
    return mDistanceMetric;
}

std::vector<double *>
DataModelingConfigurations::GetVariance(){
    return mVariance;
}

bool
DataModelingConfigurations::GetLog(){
    return mLog;
}

void
DataModelingConfigurations::SetLog(bool aLog){
     mLog = aLog;
}

void
//should datalog be smth like "verbose"??
DataModelingConfigurations::ParseDataLog(std::string aLog){
    if(aLog == "yes" || aLog == "Yes"){
        SetLog(true);
    }else if(aLog == "no" || aLog == "No"){
        SetLog(false);
    }else{
        throw range_error("Invalid value. Please use yes or no values only.");
    }
}

FILE *
DataModelingConfigurations::GetFileLog(){
    if (this->mpFileLog == nullptr) {
        throw std::runtime_error("Log File is null");
    }
    return mpFileLog;
}

double
DataModelingConfigurations::GetAvgExecutedTimePerIteration(){
    return mAvgExecutedTimePerIteration;
}

void
DataModelingConfigurations::SetAvgExecutedTimePerIteration(double aAvgExecTimePerIter){
    mAvgFlopsPerIteration = aAvgExecTimePerIter;
}

double
DataModelingConfigurations::GetAvgFlopsPerIteration(){
    return mAvgFlopsPerIteration;
}

void
DataModelingConfigurations::SetAvgFlopsPerIteration(double aAvgFlopsPerIter){
    mAvgFlopsPerIteration = aAvgFlopsPerIter;
}

double
DataModelingConfigurations::GetFinalLogLik(){
    return mFinalLogLik;
}

void
DataModelingConfigurations::SetFinalLogLik(double aAvgExecTimePerIter){
    mFinalLogLik = aAvgExecTimePerIter;
}


