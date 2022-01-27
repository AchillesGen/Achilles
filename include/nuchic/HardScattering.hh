#ifndef HARD_SCATTERING_HH
#define HARD_SCATTERING_HH

#include "nuchic/Spinor.hh"
#include <utility>
#include <complex>
#include <vector>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "yaml-cpp/yaml.h"
#pragma GCC diagnostic pop

#include "nuchic/HardScatteringEnum.hh"
#include "nuchic/Beams.hh"
#include "nuchic/RunModes.hh"
#include "nuchic/FormFactor.hh"
#include "nuchic/Histogram.hh"
#include "nuchic/ProcessInfo.hh"
#include "plugins/Sherpa/SherpaMEs.hh"

namespace nuchic {

class FourVector;
class Particle;
class Nucleus;
class Event;
class NuclearModel;

using Particles = std::vector<Particle>;
using Current = std::vector<std::vector<std::complex<double>>>;
using Currents = std::map<int, Current>;

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
        void SetProcess(const nuchic::Process_Info&);
        nuchic::Process_Info Process() const { return m_leptonicProcess; }

        // Pointer operations
        void SetSherpa(SherpaMEs *const _sherpa) { p_sherpa = _sherpa; }
        void SetNuclear(std::unique_ptr<NuclearModel> model) { m_nuclear = std::move(model); }
        NuclearModel* Nuclear() { return m_nuclear.get(); }

    private:
        Currents LeptonicCurrents(const std::vector<FourVector>&,
                                  const double&) const;

        SherpaMEs *p_sherpa{nullptr};
        nuchic::Process_Info m_leptonicProcess;
        bool m_fill{false};
        std::unique_ptr<NuclearModel> m_nuclear{};
};

}

#endif
