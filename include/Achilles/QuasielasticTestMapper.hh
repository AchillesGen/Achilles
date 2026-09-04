// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef QUASIELASTIC_TEST_MAPPER_HH
#define QUASIELASTIC_TEST_MAPPER_HH

#include "Achilles/Beams.hh"
#include "Achilles/Mapper.hh"
#include "Achilles/RunModes.hh"
#include "Achilles/UnitsIO.hh"
#include <cmath>
#include <iostream>

namespace YAML {
class Node;
}

namespace achilles {

template <int D> class FourVectorT;
using FourVector = FourVectorT<1>;
class Beam;

/// @deprecated Not referenced anywhere in the generator; kept only for its
/// tests. Slated for removal -- prefer the standard phase-space mappers.
/// Its FixedAngleEnergy mode reports d(sigma)/dE dOmega in nb/(MeV sr) by
/// normalising to a one degree window rather than adding the Jacobian to the
/// weight, so its numbers are not directly comparable to the other mappers.
class QuasielasticTestMapper : public Mapper<FourVector> {
  public:
    QuasielasticTestMapper(const YAML::Node &, std::shared_ptr<Beam>);

    void GeneratePoint(std::vector<FourVector> &, const std::vector<double> &) override;
    double GenerateWeight(const std::vector<FourVector> &, std::vector<double> &) override;
    size_t NDims() const override { return nvars; }

  private:
    RunMode mode;
    size_t nvars;
    double m_angle;
    // Carries the one degree window factor, so it is not a bare energy.
    double m_lepton_energy;
    std::shared_ptr<Beam> m_beam;
    static constexpr double dPhi = 2 * M_PI;
    static constexpr double dCos = 2;
    static constexpr double dp = 800;
};

} // namespace achilles

#endif
