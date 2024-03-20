#ifndef EVENTGEN_HH
#define EVENTGEN_HH

#include "Achilles/CombinedCuts.hh"
#include "Achilles/Histogram.hh"
#include "Achilles/ParticleInfo.hh"
#include "Achilles/QuasielasticTestMapper.hh"
#include "Achilles/Vegas.hh"
#include "Achilles/MultiChannel.hh"
#include "Achilles/Unweighter.hh"

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
        EventGen(const std::string&, std::vector<std::string>);
        void Initialize();
        void GenerateEvents();

    private:
        bool runCascade{false}, outputEvents{false}, doHardCuts{false}, doEventCuts{false};
        bool runDecays{true};
        bool doRotate{false};
        double GenerateEvent(const std::vector<FourVector>&, const double&);
        bool MakeCuts(Event&);
        // bool MakeEventCuts(Event&);
        void Rotate(Event&);

        std::shared_ptr<Beam> beam;
        std::shared_ptr<Nucleus> nucleus;
        std::shared_ptr<Cascade> cascade;
        std::shared_ptr<HardScattering> scattering;
        CutCollection hard_cuts{};
        // CutCollection event_cuts{};
        MultiChannel integrator;
        Integrand<FourVector> integrand;
        YAML::Node config;

        std::ofstream outputfile;

        std::shared_ptr<EventWriter> writer;
        StatsData polarization0, polarization1;
        SherpaInterface *p_sherpa;
        std::unique_ptr<Unweighter> unweighter;
        std::chrono::time_point<std::chrono::system_clock> now, prev;
};

}

#endif
