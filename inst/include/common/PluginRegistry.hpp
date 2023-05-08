// Copyright (c) 2017-2023 King Abdullah University of Science and Technology,
// Copyright (C) 2023 by Brightskies inc,
// All rights reserved.
// ExaGeoStat is a software package, provided by King Abdullah University of Science and Technology (KAUST).

/**
 * @file PluginRegistry.hpp
 * @brief 
 * @version 1.0.0
 * @author Sameh Abdulah
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
        template<typename T>
        class PluginRegistry {
        public:

            typedef std::function<T *()> FactoryFunction;
            typedef std::unordered_map<std::string, FactoryFunction> FactoryMap;

            static bool Add(const std::string &name, FactoryFunction fac) {
                auto map = GetFactoryMap();
                if (map.find(name) != map.end()) {
                    return false;
                }
                GetFactoryMap()[name] = fac;
                return true;
            }

            static T *Create(const std::string &name) {
                auto map = GetFactoryMap();

                if (map.find(name) == map.end()) {
                    return nullptr;
                }
                return map[name]();
            }

        private:
            // Use Meyer's singleton to prevent SIOF
            static FactoryMap &GetFactoryMap() {
                static FactoryMap map;
                return map;
            }
        };
    }//namespace common
}//namespace exageostat

#endif //EXAGEOSTATCPP_PLUGINREGISTRY_HPP
