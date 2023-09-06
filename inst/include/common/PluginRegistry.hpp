
// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file PluginRegistry.hpp
 * @brief Defines a template class for registering and creating plugins.
 * @version 1.0.0
 * @author Mahmoud ElKarargy
 * @date 2023-04-30
**/

#ifndef EXAGEOSTATCPP_PLUGINREGISTRY_HPP
#define EXAGEOSTATCPP_PLUGINREGISTRY_HPP

#include <iostream>
#include <string>
#include <map>
#include <set>
#include <functional>

namespace exageostat {
    namespace plugins {
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
            static T *Create(const std::string &aName) {
                auto map = GetFactoryMap();

                if (map.find(aName) == map.end()) {
                    return nullptr;
                }
                return map[aName]();
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
    }//namespace plugins
}//namespace exageostat

#endif //EXAGEOSTATCPP_PLUGINREGISTRY_HPP