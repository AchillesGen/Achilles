#ifndef FINTERACTIONS_HH
#define FINTERACTIONS_HH

#include "Achilles/Interactions.hh"

extern "C" {
    using namespace achilles;

    // Procedures implemented in Fortran to interface with Interactions
    void InitializeInteraction(const char*);
    double CrossSectionFortran(const Particle *part1, const Particle *part2);
    ThreeVector* MakeMomentumFortran(bool, const double, const double[2]);
}

namespace achilles {

class FortranInteraction : public Interactions {
    public:
        FortranInteraction(const YAML::Node &node) {
            auto name = node["Name"].as<std::string>();
            size_t len = name.size();
            auto cname = std::unique_ptr<char>(new char[len]);
            strcpy(cname.get(), name.c_str());
            InitializeInteraction(cname.get());
        };

        static std::unique_ptr<Interactions> Create(const YAML::Node &data) {
            return std::make_unique<FortranInteraction>(data);
        }

        static std::string GetName() { return "FortranInteraction"; }
        std::string Name() const override { return FortranInteraction::GetName(); }
        static bool IsRegistered() noexcept { return registered; }
        double CrossSection(const Particle &part1, const Particle &part2) const override {
            return CrossSectionFortran(&part1, &part2);
        }
        ThreeVector MakeMomentum(bool samePID, const double &pcm,
                                 const std::array<double, 2> &rans) const override {
            return {*MakeMomentumFortran(samePID, pcm, rans.data())};
        }

    private:
        static bool registered;
};

}

#endif
