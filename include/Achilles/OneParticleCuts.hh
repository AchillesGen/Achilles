// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef ONE_PARTICLE_CUT_HH
#define ONE_PARTICLE_CUT_HH

#include "Achilles/Cuts.hh"

namespace achilles {

/// Registry base. The bounds live in each concrete cut, because their unit
/// depends on what the cut is cutting on: an energy cut is in MeV/GeV, an
/// angular cut is a plain number of degrees.
class OneParticleCut {
  public:
    // Calculate cut for input FourVector
    virtual bool MakeCut(const FourVector &) const = 0;
    OneParticleCut() = default;
    virtual ~OneParticleCut() = default;
    static std::string Name() { return "Single Particle"; }
};

#define ONE_PARTICLE_CUT(CutName, CutType)                                         \
    class CutName##Cut : public OneParticleCut,                                    \
                         public CutBase<CutType>,                                  \
                         RegistrableCut<OneParticleCut, CutName##Cut> {            \
      public:                                                                      \
        CutName##Cut(const YAML::Node &node) : CutBase<CutType>(node) {}           \
        static std::string Name() { return #CutName; }                             \
        static std::unique_ptr<OneParticleCut> Construct(const YAML::Node &node) { \
            return std::make_unique<CutName##Cut>(node);                           \
        }                                                                          \
        bool MakeCut(const FourVector &) const override;                           \
    }

ONE_PARTICLE_CUT(Energy, units::Energy);
ONE_PARTICLE_CUT(Momentum, units::Energy);
ONE_PARTICLE_CUT(TransverseMomentum, units::Energy);
// E * theta^2 is an energy; theta is in radians.
ONE_PARTICLE_CUT(ETheta2, units::Energy);
// Angles are dimensionless; this cut is in degrees.
ONE_PARTICLE_CUT(AngleTheta, double);

// class EnergyCut : public OneParticleCut, RegistrableCut<OneParticleCut, EnergyCut> {
//     public:
//         EnergyCut(const YAML::Node &node) : OneParticleCut(node) {}
//
//         static std::string Name() { return "Energy"; }
//         static std::unique_ptr<OneParticleCut> Construct(const YAML::Node &node) {
//             return std::make_unique<EnergyCut>(node);
//         }
//
//         bool MakeCut(const FourVector&) const override { return true; }
//
//     private:
//         cut_range m_range;
// };

// class MomentumCut : public OneParticleCut, RegistrableCut<OneParticleCut, MomentumCut> {
//     public:
//         MomentumCut() = default;
//
//         static std::string Name() { return "Momentum"; }
//         static std::unique_ptr<OneParticleCut> Construct(const YAML::Node&) {
//             return std::make_unique<MomentumCut>();
//         }
//
//         bool MakeCut(const FourVector&) const override { return true; }
//
//     private:
//         cut_range m_range;
// };
//
// class AngleThetaCut : public OneParticleCut, RegistrableCut<OneParticleCut, AngleThetaCut> {
//     public:
//         AngleThetaCut() = default;
//
//         static std::string Name() { return "AngleTheta"; }
//         static std::unique_ptr<OneParticleCut> Construct(const YAML::Node&) {
//             return std::make_unique<AngleThetaCut>();
//         }
//
//         bool MakeCut(const FourVector&) const override { return true; }
//
//     private:
//         cut_range m_range;
// };
//
// class TransverseMomentumCut : public OneParticleCut, RegistrableCut<OneParticleCut,
// TransverseMomentumCut> {
//     public:
//         TransverseMomentumCut() = default;
//
//         static std::string Name() { return "TransverseMomentum"; }
//         static std::unique_ptr<OneParticleCut> Construct(const YAML::Node&) {
//             return std::make_unique<TransverseMomentumCut>();
//         }
//
//         bool MakeCut(const FourVector&) const override { return true; }
//
//     private:
//         cut_range m_range;
// };
//
// class ETheta2Cut : public OneParticleCut, RegistrableCut<OneParticleCut, ETheta2Cut> {
//     public:
//         ETheta2Cut() = default;
//
//         static std::string Name() { return "ETheta2"; }
//         static std::unique_ptr<OneParticleCut> Construct(const YAML::Node&) {
//             return std::make_unique<ETheta2Cut>();
//         }
//
//         bool MakeCut(const FourVector&) const override { return true; }
//
//     private:
//         cut_range m_range;
// };

} // namespace achilles

#endif
