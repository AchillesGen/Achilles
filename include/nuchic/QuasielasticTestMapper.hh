#ifndef QUASIELASTIC_TEST_MAPPER_HH
#define QUASIELASTIC_TEST_MAPPER_HH

#include <cmath>
#include <iostream>
#include "nuchic/Mapper.hh"
#include "nuchic/RunModes.hh"
#include "nuchic/Beams.hh"

namespace YAML {
    class Node;
}

namespace nuchic {

class FourVector;
class Beam;

class QuasielasticTestMapper : public Mapper<FourVector> {
    public:
        QuasielasticTestMapper(const YAML::Node&, std::shared_ptr<Beam>);

        void GeneratePoint(std::vector<FourVector>&, const std::vector<double>&) override;
        double GenerateWeight(const std::vector<FourVector>&, std::vector<double>&) override;
        size_t NDims() const override { return nvars; }
        YAML::Node ToYAML() const override {
            YAML::Node result;
            result["Name"] = "QuasielasticTest";
            result["Beam"] = *m_beam;
            return result;
        }

    private:
        RunMode mode;
        size_t nvars;
        double m_angle, m_lepton_energy;
        std::shared_ptr<Beam> m_beam;
        static constexpr double dPhi = 2*M_PI;
        static constexpr double dCos = 2;
        static constexpr double dp = 800;
};

}

#endif
