// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef CUTS_HH
#define CUTS_HH

#include <utility>

#include "Achilles/Factory.hh"
#include "Achilles/ParticleInfo.hh"
#include "Achilles/UnitsFormat.hh"
#include "Achilles/UnitsIO.hh"
#include "fmt/ranges.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "yaml-cpp/yaml.h"
#pragma GCC diagnostic pop

namespace achilles {

template <int D> class FourVectorT;
using FourVector = FourVectorT<1>;

/// Widest representable cut bound, for both plain numbers and Quantities.
template <class T, bool = std::is_arithmetic<T>::value> struct CutLimits;
template <class T> struct CutLimits<T, true> {
    static T lowest() { return std::numeric_limits<T>::lowest(); }
    static T max() { return std::numeric_limits<T>::max(); }
};
template <class T> struct CutLimits<T, false> {
    static T lowest() { return T{std::numeric_limits<double>::lowest()}; }
    static T max() { return T{std::numeric_limits<double>::max()}; }
};

template <class T> class CutBase {
  public:
    using cut_range = std::pair<T, T>;
    using cut_ranges = std::vector<std::pair<T, T>>;
    CutBase() = default;

    CutBase(const YAML::Node &node) {
        if(node["min"] && node["max"] && node["range"]) {
            throw std::runtime_error("CutRange: Invalid syntax. Got min, max, and range");
        } else if(node["min"] || node["max"]) {
            T min = CutLimits<T>::lowest();
            T max = CutLimits<T>::max();
            // A dimensionful bound is a { value:, unit: } map, so only a plain
            // number is required to be a scalar here.
            if(node["min"]) {
                if(std::is_arithmetic<T>::value && !node["min"].IsScalar())
                    throw std::runtime_error("CutRange: Invalid min value");
                min = node["min"].as<T>();
            }
            if(node["max"]) {
                if(std::is_arithmetic<T>::value && !node["max"].IsScalar())
                    throw std::runtime_error("CutRange: Invalid max value");
                max = node["max"].as<T>();
            }
            m_range = {{min, max}};
        } else if(node["range"]) {
            m_range = node["range"].as<cut_ranges>();
        } else {
            throw std::runtime_error("CutRange: Invalid syntax. Missing cut values");
        }
        if(m_range.size() == 1) {
            spdlog::trace("Found cut range: [{}, {}]", m_range[0].first, m_range[0].second);
        } else {
            std::string ranges{};
            for(const auto &range : m_range) {
                ranges += fmt::format("[{}, {}], ", range.first, range.second);
            }
            spdlog::trace("Found cut range: [{}]", ranges.substr(0, ranges.size() - 2));
        }
    }

  protected:
    bool CheckCut(const T &val) const {
        bool result = false;
        // Run through possible cut ranges, keeping the event if a single cut range
        // is satisfied. Note that these cut ranges are all for a single variable,
        // e.g., theta in [0,1] OR [2,3] OR ...
        // Return true if any single cut range is satisfied.
        for(const auto &cut : m_range) { result |= (val > cut.first && val < cut.second); }

        return result;
    }

  private:
    cut_ranges m_range;
};

template <class Base, class Derived>
using RegistrableCut = Registrable<Base, Derived, const YAML::Node &>;
template <class Base> using CutFactory = Factory<Base, const YAML::Node &>;

} // namespace achilles

namespace YAML {

// A bound is either a plain number or a { value:, unit: } map, so these check
// the shape of the sequence and leave the element conversion to convert<T>.
template <typename T> struct convert<std::pair<T, T>> {
    static bool decode(const Node &node, std::pair<T, T> &range) {
        if(node.IsSequence() && node.size() == 2) {
            range = {node[0].as<T>(), node[1].as<T>()};
            return true;
        }
        return false;
    }
};

template <typename T> struct convert<std::vector<std::pair<T, T>>> {
    static bool decode(const Node &node, std::vector<std::pair<T, T>> &cutRange) {
        if(!node.IsSequence()) return false;
        if(node.size() == 2 && !node[0].IsSequence()) {
            cutRange = {{node[0].as<T>(), node[1].as<T>()}};
        } else {
            for(const auto &subNode : node) { cutRange.push_back(subNode.as<std::pair<T, T>>()); }
        }
        return true;
    }
};

} // namespace YAML

#endif
