
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file PluginRegistry.hpp
 * @brief Defines a template class for registering and creating plugins.
 * @version 1.0.1
 * @author Mahmoud ElKarargy
 * @date 2023-04-30
**/

#ifndef EXAGEOSTATCPP_PLUGINREGISTRY_HPP
#define EXAGEOSTATCPP_PLUGINREGISTRY_HPP

//// TODO: This is a hot fix to avoid the problem in HiCMA which set the min definition with a conflict implementation of chrono library.
#ifdef min
#undef min
#endif

#include <iostream>
#include <string>
#include <map>
#include <set>
#include <functional>

#include <configurations/Configurations.hpp>

namespace exageostat::plugins {
    /**
     * @brief Template class for registering and creating plugins.
     * @tparam T Data Type: float or double
     *
     */
    template<typename T>
    class PluginRegistry {
    public:
        /**
         * @brief Function type that returns a pointer to an instance of T.
         *
         */
        typedef std::function<T *()> FactoryFunction;

        /**
         * @brief Unordered map that maps plugin names to their corresponding factory functions.
         *
         */
        typedef std::unordered_map<std::string, FactoryFunction> FactoryMap;

        /**
         * @brief Adds a factory function to the FactoryMap under the given plugin name.
         * @param[in] name The name of the plugin to be added.
         * @param[in] fac The factory function to be added.
         * @return true if the factory function was successfully added, false otherwise.
         *
         */
        static bool Add(const std::string &name, FactoryFunction fac) {
            auto map = GetFactoryMap();
            if (map.find(name) != map.end()) {
                return false;
            }
            GetFactoryMap()[name] = fac;
            return true;
        }

        /**
         * @brief Creates an instance of the plugin with the given name.
         * @param[in] name The name of the plugin to be created.
         * @return A pointer to the created plugin, or nullptr if the plugin could not be created.
         *
         */
        static T *Create(const std::string &aName, const int& aTimeSlot) {
            auto map = GetFactoryMap();

            if (map.find(aName) == map.end()) {
                return nullptr;
            }
            // Get the object from the map.
            T* object = map[aName]();
            // Automatically set the new P value which will get updated with the user input of P.
            object->SetPValue(aTimeSlot);
            return object;
        }

    private:
        /**
         * @brief Returns a reference to the FactoryMap singleton.
         * @return A reference to the FactoryMap singleton.
         *
         */
        static FactoryMap &GetFactoryMap() {
            static FactoryMap mSelfRegisteringMap;
            return mSelfRegisteringMap;
        }
    };
}//namespace exageostat

#endif //EXAGEOSTATCPP_PLUGINREGISTRY_HPP