#ifndef FINTERACTIONS_HH
#define FINTERACTIONS_HH

#include "nuchic/Interactions.hh"

extern "C" {
    using namespace nuchic;

    // Procedures implemented in Fortran to interface with Interactions
    void InitializeInteraction(const char*);
    double CrossSectionFortran(const Particle *part1, const Particle *part2);
    ThreeVector* MakeMomentumFortran(bool, const double, const double[2]);
}

namespace nuchic {

class FortranInteraction : public Interactions {
    public:
        FortranInteraction(const std::string &name) {
            size_t len = name.size();
            char* cname = new char[len];
            strcpy(cname, name.c_str());
            fmt::print("{}, {}\n", cname, strlen(cname));
            InitializeInteraction(cname);
            delete[] cname;
        };

        static std::unique_ptr<Interactions> Create(const std::string &data) {
            return std::make_unique<FortranInteraction>(data);
        }

        static std::string GetName() { return "FortranInteraction"; }
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
