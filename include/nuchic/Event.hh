#ifndef EVENT_HH
#define EVENT_HH

#include <map>
#include <memory>
#include <utility>
#include <vector>

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

class Event {
    struct PhaseSpaceStruct {
        std::vector<FourVector> momentum;
        double weight{};
    };

    struct MatrixElementStruct {
        std::vector<PID> inital_state;
        std::vector<PID> final_state;
        double weight{};
    };

    using MatrixElementVec = std::vector<MatrixElementStruct>;

    public:
        Event() = default;
        Event(std::shared_ptr<Nucleus>, std::shared_ptr<Beam>,
              const std::vector<double>&, double);

        void SetHardScatteringType(HardScatteringType type) { m_type = type; }
        void InitializeLeptons(size_t);
        void InitializeHadrons(const std::vector<std::array<size_t, 3>>&);
        void Finalize();

        const NuclearRemnant &Remnant() const { return m_remnant; }

        const PhaseSpaceStruct &PhaseSpace() const { return m_ps; }
        PhaseSpaceStruct &PhaseSpace() { return m_ps; }

        const MatrixElementStruct &MatrixElement(size_t i) const { return m_me[i]; }
        MatrixElementStruct &MatrixElement(size_t i) { return m_me[i]; }

        const MatrixElementVec &MatrixElements() const { return m_me; }
        MatrixElementVec &MatrixElements() { return m_me; }

        bool TotalCrossSection();
        std::vector<double> EventProbs() const;
        static bool MatrixCompare(const MatrixElementStruct&, double);

        const std::shared_ptr<Nucleus>& CurrentNucleus() const { return m_nuc; }
        std::shared_ptr<Nucleus>& CurrentNucleus() { return m_nuc; }
        const std::shared_ptr<Beam>& CurrentBeam() const { return m_beam; }
        std::shared_ptr<Beam>& CurrentBeam() { return m_beam; }

        void AddParticle(const Particle&);
        vParticles Particles() const;
        const vParticles& Hadrons() const;
        vParticles& Hadrons();
        const vParticles& Leptons() const { return m_leptons; }
        vParticles& Leptons() { return m_leptons; }
        double Weight() const;

    private:
        static double AddEvents(double, const MatrixElementStruct&);

        bool ValidateEvent(size_t) const;

        HardScatteringType m_type{HardScatteringType::None};
        std::shared_ptr<Nucleus> m_nuc;
        NuclearRemnant m_remnant{};
        std::shared_ptr<Beam> m_beam;
        PhaseSpaceStruct m_ps{};
        MatrixElementVec m_me;
        double m_vWgt{}, m_meWgt{};
        vParticles m_leptons{};
};

}

#endif
