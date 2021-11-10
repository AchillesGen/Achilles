#ifndef EVENT_HH
#define EVENT_HH

#include <map>
#include <memory>
#include <utility>
#include <vector>

#include "nuchic/Nuchic.hh"
#include "nuchic/HardScatteringEnum.hh"
#include "nuchic/NuclearRemnant.hh"

namespace nuchic {

class PID;
class FourVector;
class Particle;
class Nucleus;
class Beam;

struct InitialState {
    std::shared_ptr<Beam> beam;
    std::shared_ptr<Nucleus> nucleus;
};

using vParticles = std::vector<Particle>;
using vMomentum = std::vector<FourVector>;

class Event {
    struct MatrixElementStruct {
        std::vector<PID> inital_state;
        std::vector<PID> final_state;
        double weight{};
    };

    using MatrixElementVec = std::vector<MatrixElementStruct>;

    public:
        Event() = default;
        Event(std::shared_ptr<Nucleus>, 
              std::vector<FourVector>, double);
        MOCK ~Event() = default;

        void SetHardScatteringType(HardScatteringType type) { m_type = type; }
        void InitializeLeptons(size_t);
        void InitializeHadrons(const std::vector<std::array<size_t, 3>>&);
        void InitializeCoherent();
        void Finalize();

        const NuclearRemnant &Remnant() const { return m_remnant; }

        const vMomentum &Momentum() const { return m_mom; }
        vMomentum &Momentum() { return m_mom; }

        const MatrixElementStruct &MatrixElement(size_t i) const { return m_me[i]; }
        MatrixElementStruct &MatrixElement(size_t i) { return m_me[i]; }

        const MatrixElementVec &MatrixElements() const { return m_me; }
        MatrixElementVec &MatrixElements() { return m_me; }

        bool TotalCrossSection();
        std::vector<double> EventProbs() const;
        static bool MatrixCompare(const MatrixElementStruct&, double);

        const std::shared_ptr<Nucleus>& CurrentNucleus() const { return m_nuc; }
        MOCK std::shared_ptr<Nucleus>& CurrentNucleus() { return m_nuc; }

        void AddParticle(const Particle&);
        vParticles Particles() const;
        const vParticles& Hadrons() const;
        MOCK vParticles& Hadrons();
        const vParticles& Leptons() const { return m_leptons; }
        vParticles& Leptons() { return m_leptons; }
        double Weight() const;
        void SetMEWeight(double wgt) { m_meWgt = wgt; }
        void Rotate(const std::array<double,9>&);

        bool IsCoherent() const { return m_coh; }
        double &CoherentXsec() { return m_xsec_coherent; }
        const double &CoherentXsec() const { return m_xsec_coherent; }

    private:
        static double AddEvents(double, const MatrixElementStruct&);

        bool ValidateEvent(size_t) const;

        HardScatteringType m_type{HardScatteringType::None};
        std::shared_ptr<Nucleus> m_nuc;
        NuclearRemnant m_remnant{};
        vMomentum m_mom{};
        MatrixElementVec m_me;
        double m_vWgt{}, m_meWgt{};
        vParticles m_leptons{};
        vParticles m_history{};

        // Coherent variables
        bool m_coh{false};
        double m_xsec_coherent{};
};

}

#endif
