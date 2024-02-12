
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file SyntheticGenerator.cpp
 * @brief Implementation of the SyntheticGenerator class
 * @version 1.0.1
 * @author Mahmoud ElKarargy
 * @author Sameh Abdulah
 * @date 2023-02-14
**/

#include <data-generators/concrete/SyntheticGenerator.hpp>
#include <data-generators/LocationGenerator.hpp>

using namespace exageostat::generators::synthetic;
using namespace exageostat::dataunits;
using namespace exageostat::common;
using namespace exageostat::configurations;

template<typename T>
SyntheticGenerator<T> *SyntheticGenerator<T>::GetInstance() {

    if (mpInstance == nullptr) {
        mpInstance = new SyntheticGenerator<T>();
    }
    return mpInstance;
}

template<typename T>
std::unique_ptr<ExaGeoStatData<T>>
SyntheticGenerator<T>::CreateData(exageostat::configurations::Configurations &aConfigurations,
                                  const exageostat::hardware::ExaGeoStatHardware &aHardware,
                                  exageostat::kernels::Kernel<T> &aKernel) {

    int n = aConfigurations.GetProblemSize() * aConfigurations.GetTimeSlot();
    auto data = std::make_unique<ExaGeoStatData<T>>(n, aConfigurations.GetDimension());
    // Allocated new Locations object.
    auto *locations = new Locations<T>(n, aConfigurations.GetDimension());

    int parameters_number = aKernel.GetParametersNumbers();

    // Set initial theta values.
    Configurations::InitTheta(aConfigurations.GetInitialTheta(), parameters_number);
    aConfigurations.SetInitialTheta(aConfigurations.GetInitialTheta());

    // Generate Locations phase
    LocationGenerator<T>::GenerateLocations(n, aConfigurations.GetTimeSlot(), aConfigurations.GetDimension(), *locations);
    data->SetLocations(*locations);

    // Generate Descriptors phase
    auto linear_algebra_solver = linearAlgebra::LinearAlgebraFactory<T>::CreateLinearAlgebraSolver(EXACT_DENSE);
    linear_algebra_solver->GenerateSyntheticData(aConfigurations, aHardware, data, aKernel);

    if (aConfigurations.GetLogger()) {
        VERBOSE("Writing generated data to the disk (Synthetic Dataset Generation Phase) .....")
#ifdef CHAMELEON_USE_MPI
        auto pMatrix = new T[aConfigurations.GetProblemSize()];
        std::string path = aConfigurations.GetLoggerPath();

        CHAMELEON_Desc2Lap(ChamUpperLower, data->GetDescriptorData()->GetDescriptor(CHAMELEON_DESCRIPTOR,
                                                                               DESCRIPTOR_Z).chameleon_desc,  pMatrix, aConfigurations.GetProblemSize());
        if (helpers::CommunicatorMPI::GetInstance()->GetRank() == 0){
            helpers::DiskWriter<T>::WriteVectorsToDisk(*((T *) data->GetDescriptorData()->GetDescriptor(CHAMELEON_DESCRIPTOR,
                                                                                               DESCRIPTOR_Z).chameleon_desc->mat),
                                              aConfigurations.GetProblemSize(), parameters_number, path,
                                              *data->GetLocations());
        }
        delete[] pMatrix;
#else
        std::string path = aConfigurations.GetLoggerPath();
        helpers::DiskWriter<T>::WriteVectorsToDisk(*((T *) data->GetDescriptorData()->GetDescriptor(CHAMELEON_DESCRIPTOR,
                                                                                           DESCRIPTOR_Z).chameleon_desc->mat),
                                          aConfigurations.GetProblemSize(), aKernel.GetVariablesNumber(), path,
                                          *data->GetLocations());
#endif
        VERBOSE("Done.")
    }
    results::Results::GetInstance()->SetGeneratedLocationsNumber(n);
    results::Results::GetInstance()->SetIsLogger(aConfigurations.GetLogger());
    results::Results::GetInstance()->SetLoggerPath(aConfigurations.GetLoggerPath());

    return data;
}

template<typename T>
void SyntheticGenerator<T>::ReleaseInstance() {
    if (mpInstance != nullptr) {
        mpInstance = nullptr;
    }
}

template<typename T> SyntheticGenerator<T> *SyntheticGenerator<T>::mpInstance = nullptr;
