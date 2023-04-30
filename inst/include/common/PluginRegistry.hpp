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
#include <kernels/Kernel.hpp>

#define EXAGEOSTAT_REGISTER_PLUGIN(plugin_name, create_func) \
    bool plugin_name ## _entry = plugins::PluginRegistry<exageostat::kernels::Kernel>::add(#plugin_name, (create_func))


namespace exageostat {
    namespace plugins {
        template<typename T>
        class PluginRegistry {
        public:

            typedef std::function<T *()> FactoryFunction;
            typedef std::unordered_map <std::string, FactoryFunction> FactoryMap;

            static bool add(const std::string &name, FactoryFunction fac) {
                auto map = getFactoryMap();
                std::cout << "bt add bs?" << std::endl;

                if (map.find(name) != map.end()) {
                    return false;
                }

                std::cout << "bt add bs?" << name << std::endl;
                getFactoryMap()[name] = fac;
                return true;
            }

            static T *create(const std::string &name) {
                auto map = getFactoryMap();
                std::cout << "malk bs?" << std::endl;
                if (map.find(name) == map.end()) {
                    std::cout << "LEEEEEEEEEEEH" << std::endl;
                    return nullptr;
                }

                return map[name]();
            }

        private:
            // Use Meyer's singleton to prevent SIOF
            static FactoryMap &getFactoryMap() {
                static FactoryMap map;
                return map;
            }
        };
    }//namespace common
}//namespace exageostat

#endif //EXAGEOSTATCPP_PLUGINREGISTRY_HPP
