#ifndef FINTERACTIONS_HH
#define FINTERACTIONS_HH

#include "Achilles/Interactions.hh"
#include "Achilles/Particle.hh"

extern "C" {
using namespace achilles;

// Procedures implemented in Fortran to interface with Interactions
void InitializeInteraction(const char *);
double CrossSectionFortran(const Particle *part1, const Particle *part2);
ThreeVector *MakeMomentumFortran(bool, const double, const double[2]);
}

namespace achilles {

class FortranInteraction : public Interaction, RegistrableInteraction<FortranInteraction> {
  public:
    FortranInteraction(const YAML::Node &node) {
        throw std::runtime_error("FortranInteraction is not implemented");
        auto name = node["Name"].as<std::string>();
        size_t len = name.size();
        auto cname = std::unique_ptr<char>(new char[len]);
        strcpy(cname.get(), name.c_str());
        InitializeInteraction(cname.get());
    };

    static std::unique_ptr<Interaction> Construct(const YAML::Node &data) {
        return std::make_unique<FortranInteraction>(data);
    }

    static std::string Name() { return "FortranInteraction"; }
    std::string GetName() const override { return FortranInteraction::GetName(); }
    static bool IsRegistered() noexcept { return registered; }
    InteractionResults CrossSection(const Particle &part1, const Particle &part2) const override {
        return {{{part1.ID(), part2.ID()}, CrossSectionFortran(&part1, &part2)}};
    }
    std::vector<std::pair<int, int>> InitialStates() const override {
        return {{2112, 2112}, {2112, 2212}, {2212, 2212}};
    }
    std::vector<Particle> GenerateMomentum(const Particle &, const Particle &,
                                           const std::vector<PID> &,
                                           const std::vector<double> &) const override {
        std::vector<Particle> particles;
        return particles;
    }
    ThreeVector MakeMomentum(bool samePID, const double &pcm,
                             const std::array<double, 2> &rans) const {
        return {*MakeMomentumFortran(samePID, pcm, rans.data())};
    }

  private:
    static bool registered;
};

} // namespace achilles

#endif
