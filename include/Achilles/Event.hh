#ifndef EVENT_HH
#define EVENT_HH

#include <map>
#include <memory>
#include <utility>
#include <vector>

#include "Achilles/Achilles.hh"
#include "Achilles/HardScatteringEnum.hh"
#include "Achilles/NuclearRemnant.hh"
#include "Achilles/ProcessInfo.hh"
#include "Achilles/EventHistory.hh"

namespace achilles {

class PID;
class FourVector;
class Particle;
class Nucleus;
class Beam;
class NuclearModel;

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
        Event(double vWgt = 0) : m_vWgt{vWgt} {}
        Event(std::shared_ptr<Nucleus>, 
              std::vector<FourVector>, double);
        MOCK ~Event() = default;

        void SetHardScatteringType(HardScatteringType type) { m_type = type; }
        MOCK void InitializeLeptons(const Process_Info&);
        MOCK void InitializeHadrons(const Process_Info&);
        void Finalize();

        MOCK const NuclearRemnant &Remnant() const { return m_remnant; }

        MOCK const vMomentum &Momentum() const { return m_mom; }
        MOCK vMomentum &Momentum() { return m_mom; }

        MOCK const EventHistory &History() const { return m_history; }
        EventHistory &History() { return m_history; }

        const double &MatrixElementWgt(size_t i) const { return m_me[i]; }
        double &MatrixElementWgt(size_t i) { return m_me[i]; }

        const std::vector<double> &MatrixElementWgts() const { return m_me; }
        std::vector<double> &MatrixElementWgts() { return m_me; }

        bool TotalCrossSection();
        size_t SelectNucleon() const;

        MOCK const std::shared_ptr<Nucleus> CurrentNucleus() const { return m_nuc; }
        MOCK std::shared_ptr<Nucleus> CurrentNucleus() { return m_nuc; }

        const double& Flux() const { return flux; }
        double& Flux() { return flux; }

        MOCK vParticles Particles() const;
        MOCK const vParticles& Hadrons() const;
        MOCK vParticles& Hadrons();
        MOCK const vParticles& Leptons() const { return m_leptons; }
        MOCK vParticles& Leptons() { return m_leptons; }
        void CalcWeight();
        MOCK const double& Weight() const { return m_wgt; }
        MOCK double& Weight() { return m_wgt; }
        void SetMEWeight(double wgt) { m_meWgt = wgt; }
        void Rotate(const std::array<double,9>&);

        bool operator==(const Event &other) const {
            return m_type == other.m_type && m_nuc == other.m_nuc
                && m_remnant == other.m_remnant && m_mom == other.m_mom
                && m_me == other.m_me && m_vWgt == other.m_meWgt
                && m_leptons == other.m_leptons;
        }

    private:
        static bool MatrixCompare(const MatrixElementStruct&, double);
        static double AddEvents(double, const MatrixElementStruct&);
        std::vector<double> EventProbs() const;

        // bool ValidateEvent(size_t) const;

        HardScatteringType m_type{HardScatteringType::None};
        std::shared_ptr<Nucleus> m_nuc;
        NuclearRemnant m_remnant{};
        vMomentum m_mom{};
        std::vector<double> m_me;
        double m_vWgt{}, m_meWgt{}, m_wgt{-1};
        vParticles m_leptons{};
        EventHistory m_history{};
        double flux;
};

}

#endif
