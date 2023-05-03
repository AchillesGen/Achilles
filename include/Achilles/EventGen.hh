#ifndef EVENTGEN_HH
#define EVENTGEN_HH

#include "Achilles/CombinedCuts.hh"
#include "Achilles/Histogram.hh"
#include "Achilles/MultiChannel.hh"
#include "Achilles/ParticleInfo.hh"
#include "Achilles/Process.hh"
#include "Achilles/QuasielasticTestMapper.hh"
#include "Achilles/Unweighter.hh"
#include "Achilles/Vegas.hh"

#include <memory>
#include <vector>

namespace YAML {
class Node;
}

namespace achilles {

// Forward declare types
class Event;
class Beam;
class Nucleus;
class Cascade;
class HardScattering;
class EventWriter;

class SherpaInterface;

class EventGen {
  public:
    EventGen(const std::string &, std::vector<std::string>);
    void Initialize();
    void GenerateEvents();

  private:
    bool GenerateSingleEvent();

    bool runCascade{false}, outputEvents{false}, doHardCuts{false}, doEventCuts{false};
    bool runDecays{false};
    bool doRotate{false};
    double GenerateEvent(const std::vector<FourVector> &, const double &);
    bool MakeCuts(Event &);
    // bool MakeEventCuts(Event&);
    void Rotate(Event &);

    std::shared_ptr<Beam> beam;
    std::shared_ptr<Nucleus> nucleus;
    std::shared_ptr<Cascade> cascade;
    std::shared_ptr<HardScattering> scattering;
    std::vector<ProcessGroup> process_groups;
    CutCollection hard_cuts{};
    // CutCollection event_cuts{};
    MultiChannel integrator;
    Integrand<FourVector> integrand;
    YAML::Node config;
    std::vector<double> m_group_weights{};
    double m_max_weight{};

    std::shared_ptr<EventWriter> writer;
    std::unique_ptr<Unweighter> unweighter;
    SherpaInterface *p_sherpa;
};

} // namespace achilles

#endif
