#ifndef QUASIELASTIC_TEST_MAPPER_HH
#define QUASIELASTIC_TEST_MAPPER_HH

#include "Achilles/Beams.hh"
#include "Achilles/Mapper.hh"
#include "Achilles/RunModes.hh"
#include <cmath>
#include <iostream>

namespace YAML {
class Node;
}

namespace achilles {

class FourVector;
class Beam;

class QuasielasticTestMapper : public Mapper<FourVector> {
  public:
    QuasielasticTestMapper(const YAML::Node &, std::shared_ptr<Beam>);

    void GeneratePoint(std::vector<FourVector> &, const std::vector<double> &) override;
    double GenerateWeight(const std::vector<FourVector> &, std::vector<double> &) override;
    size_t NDims() const override { return nvars; }

  private:
    RunMode mode;
    size_t nvars;
    double m_angle, m_lepton_energy;
    std::shared_ptr<Beam> m_beam;
    static constexpr double dPhi = 2 * M_PI;
    static constexpr double dCos = 2;
    static constexpr double dp = 800;
};

} // namespace achilles

#endif
