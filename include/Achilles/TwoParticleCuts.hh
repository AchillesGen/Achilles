// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef TWO_PARTICLE_CUT_HH
#define TWO_PARTICLE_CUT_HH

#include "Achilles/Cuts.hh"

namespace achilles {

/// See OneParticleCut: the bounds live in each concrete cut so that each can
/// carry the unit of the variable it cuts on.
class TwoParticleCut {
  public:
    // Calculate cut for two input FourVectors
    virtual bool MakeCut(const FourVector &, const FourVector &) const = 0;
    TwoParticleCut() = default;
    virtual ~TwoParticleCut() = default;
    static std::string Name() { return "Two Particle"; }
};

#define TWO_PARTICLE_CUT(CutName, CutType)                                         \
    class CutName##Cut : public TwoParticleCut,                                    \
                         public CutBase<CutType>,                                  \
                         RegistrableCut<TwoParticleCut, CutName##Cut> {            \
      public:                                                                      \
        CutName##Cut(const YAML::Node &node) : CutBase<CutType>(node) {}           \
        static std::string Name() { return #CutName; }                             \
        static std::unique_ptr<TwoParticleCut> Construct(const YAML::Node &node) { \
            return std::make_unique<CutName##Cut>(node);                           \
        }                                                                          \
        bool MakeCut(const FourVector &, const FourVector &) const override;       \
    }

// Angles are dimensionless; this cut is in degrees.
TWO_PARTICLE_CUT(DeltaTheta, double);
TWO_PARTICLE_CUT(InvariantMass, units::Energy);

} // namespace achilles

#endif
