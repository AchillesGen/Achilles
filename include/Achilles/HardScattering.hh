#ifndef HARD_SCATTERING_HH
#define HARD_SCATTERING_HH

#include "Achilles/Spinor.hh"
#include <utility>
#include <complex>
#include <vector>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "yaml-cpp/yaml.h"
#pragma GCC diagnostic pop

#include "Achilles/HardScatteringEnum.hh"
#include "Achilles/Beams.hh"
#include "Achilles/RunModes.hh"
#include "Achilles/FormFactor.hh"
#include "Achilles/Histogram.hh"
#include "Achilles/ProcessInfo.hh"
#ifdef ACHILLES_SHERPA_INTERFACE
#include "plugins/Sherpa/SherpaInterface.hh"
#endif

namespace achilles {

class FourVector;
class Particle;
class Nucleus;
class Event;
class NuclearModel;

using Particles = std::vector<Particle>;
using Current = std::vector<std::vector<std::complex<double>>>;
using Currents = std::map<int, Current>;
using FFDictionary = std::map<std::pair<PID, PID>, std::vector<FormFactorInfo>>;

class LeptonicCurrent {
    public:
        LeptonicCurrent() = default;
        void Initialize(const Process_Info&);
        FFDictionary GetFormFactor();
        Currents CalcCurrents(const std::vector<FourVector>&, const double&) const;

    private:
        bool NeutralCurrent(PID, PID) const;
        bool ChargedCurrent(bool, PID, PID) const;
        std::complex<double> coupl_left{}, coupl_right{};
        double mass{}, width{};
        int pid{};
        bool anti{};
};

class HardScattering {
    public:
        HardScattering() = default;
        HardScattering(const HardScattering&) = delete;
        HardScattering(HardScattering&&) = default;
        HardScattering& operator=(const HardScattering&) = delete;
        HardScattering& operator=(HardScattering&&) = default;
        ~HardScattering() = default;

        // Calculation details
        std::vector<double> CrossSection(Event&) const;

        // Select initial state
        bool FillEvent(Event&, const std::vector<double>&) const;
        size_t SelectMatrixElement(Event&) const;

        // Process accessors
        void SetProcess(const Process_Info&);
        Process_Info Process() const { return m_leptonicProcess; }

        // Pointer operations
#ifdef ACHILLES_SHERPA_INTERFACE
        void SetSherpa(SherpaInterface *const _sherpa) { p_sherpa = _sherpa; }
#endif
        void SetNuclear(std::unique_ptr<NuclearModel> model) { m_nuclear = std::move(model); }
        NuclearModel* Nuclear() { return m_nuclear.get(); }

    private:
        Currents LeptonicCurrents(const std::vector<FourVector>&,
                                  const double&) const;
        FFDictionary SMFormFactor;

#ifdef ACHILLES_SHERPA_INTERFACE
        SherpaInterface *p_sherpa{nullptr};
#endif
        LeptonicCurrent m_current{};
        Process_Info m_leptonicProcess;
        bool m_fill{false};
        std::unique_ptr<NuclearModel> m_nuclear{};
};

}

#endif
