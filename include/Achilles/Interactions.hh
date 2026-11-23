// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include "Achilles/Factory.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/Particle.hh"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-conversion"
#include "yaml-cpp/yaml.h"
#pragma GCC diagnostic pop

namespace achilles {

class Event;
class PID;
class Particle;
class Random;

struct pid_compare {
    bool operator()(const std::pair<PID, PID> &lhs, const std::pair<PID, PID> &rhs) const;
};

struct InteractionResult {
    std::vector<PID> particles;
    double cross_section;
};

using InteractionResults = std::vector<InteractionResult>;

/// Base class for implementing interaction models. These interaction models focus on the
/// interactions that occur between hadrons during the intranuclear cascade. This base class
/// is **not** to be used to inherit from for modeling the hard interactions between leptons
/// and the nucleus.
class Interaction {
  public:
    /// @name Constructors and Destructors
    ///@{

    /// Base class constructor
    Interaction() = default;
    Interaction(const Interaction &) = default;
    Interaction(Interaction &&) = default;
    Interaction &operator=(const Interaction &) = default;
    Interaction &operator=(Interaction &&) = default;

    /// Base class constructor for classes that need data files
    // Interactions(const std::string& data) {}

    /// Default destructor
    virtual ~Interaction() = default;
    ///@}

    /// Function to list all implemented interactions between two particles
    ///@return std::vector<std::pair<PID, PID>>: List of all interactions
    virtual std::vector<std::pair<PID, PID>> InitialStates() const = 0;

    /// Function that returns the total cross-section between two particles
    ///@param kicked: The propagating particle
    ///@param remnant: The composite particle containing all background nucleons
    ///@param rem_idx: The position within the remnant of the "main" particle being interacted with
	///@param nuc: The Nucleus object within which this interaction is taking place
    ///@return double: The cross-section
    virtual double TotalCrossSection(const Particle&,Particle*,size_t,const Nucleus*) const;

    /// Function to determine the cross-section between two particles
    /// broken down into all possible final state channels
    ///@param kicked: The propagating particle
    ///@param remnant: The composite particle containing all background nucleons
    ///@param rem_idx: The position within the remnant of the "main" particle being interacted with
	///@param nuc: The Nucleus object within which this interaction is taking place
    ///@return InteractionResults: The cross-section
    virtual InteractionResults CrossSection(const Particle&,Particle*,size_t,const Nucleus*) const = 0;

    /// Function to generate the final state particles
    ///@param part1: The first particle involved with the interaction
    ///@param part2: The second particle involved with the interaction
    ///@param out_pids: The final state particles to generate
    ///@param rands: The random numbers to use in the generation
    ///@return std::vector<Particle>: The final state particles
    virtual vParticles GenerateMomentum(const Particle &part1, const Particle &part2,
                                                   const std::vector<PID> &out_pids,
                                                   Random &) const = 0;

    virtual std::string GetName() const = 0;
    static std::string Name() { return "Interaction"; }

  protected:
    double CrossSectionLab(bool, const double &) const noexcept;
};

template <typename Derived>
using RegistrableInteraction = Registrable<Interaction, Derived, const YAML::Node &>;
using InteractionFactory = Factory<Interaction, const YAML::Node &>;

} // namespace achilles

namespace YAML {
template <> struct convert<achilles::InteractionResult> {
    static bool decode(const Node &node, achilles::InteractionResult &result) {
        result.particles = node["Outgoing"].as<std::vector<achilles::PID>>();
        result.cross_section = node["CrossSection"].as<double>();
        return true;
    }
};

template <> struct convert<std::vector<achilles::InteractionResult>> {
    static bool decode(const Node &node, std::vector<achilles::InteractionResult> &results) {
        for(auto result : node) results.push_back(result.as<achilles::InteractionResult>());
        return true;
    }
};

} // namespace YAML