#ifndef HARD_SCATTERING_HH
#define HARD_SCATTERING_HH

#include "Achilles/Spinor.hh"
#include <complex>
#include <utility>
#include <vector>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "yaml-cpp/yaml.h"
#pragma GCC diagnostic pop

#include "Achilles/Beams.hh"
#include "Achilles/FormFactor.hh"
#include "Achilles/HardScatteringEnum.hh"
#include "Achilles/Histogram.hh"
#include "Achilles/LeptonicCurrent.hh"
#include "Achilles/ProcessInfo.hh"
#include "Achilles/RunModes.hh"
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
using Currents = std::map<int, Current>;
using FFDictionary = std::map<std::pair<PID, PID>, std::vector<FormFactorInfo>>;

class HardScattering {
  public:
    HardScattering() = default;
    HardScattering(const HardScattering &) = delete;
    HardScattering(HardScattering &&) = default;
    HardScattering &operator=(const HardScattering &) = delete;
    HardScattering &operator=(HardScattering &&) = default;
    ~HardScattering() = default;

    // Calculation details
    std::vector<double> CrossSection(Event &) const;

    // Select initial state
    bool FillEvent(Event &, const std::vector<double> &) const;
    size_t SelectMatrixElement(Event &) const;

    // Process accessors
    void SetProcess(const ProcessInfo &);
    ProcessInfo Process() const { return m_leptonicProcess; }

    // Pointer operations
#ifdef ACHILLES_SHERPA_INTERFACE
    void SetSherpa(SherpaInterface *const _sherpa) { p_sherpa = _sherpa; }
#endif
    void SetNuclear(std::unique_ptr<NuclearModel> model) { m_nuclear = std::move(model); }
    NuclearModel *Nuclear() { return m_nuclear.get(); }

  private:
    FFDictionary SMFormFactor;

#ifdef ACHILLES_SHERPA_INTERFACE
    SherpaInterface *p_sherpa{nullptr};
#endif
    LeptonicCurrent m_current{};
    ProcessInfo m_leptonicProcess;
    bool m_fill{false};
    std::unique_ptr<NuclearModel> m_nuclear{};
};

} // namespace achilles

#endif
