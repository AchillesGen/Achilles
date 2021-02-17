#ifndef EVENTGEN_HH
#define EVENTGEN_HH

#include "nuchic/Cuts.hh"
#include "nuchic/Histogram.hh"
#include "nuchic/ParticleInfo.hh"
#include "nuchic/Vegas.hh"

#include <memory>
#include <vector>

namespace YAML {

class Node;

}

namespace nuchic {

// Forward declare types
class Event;
class Beam;
class Nucleus;
class Cascade;
class HardScattering;
class EventWriter;

using Cuts = std::map<PID, Cut>;

class EventGen {
    public:
        EventGen(const std::string&);
        void Initialize();
        void GenerateEvents();
        size_t nevents;
        size_t total_events;
        
    private:
        bool runCascade, outputEvents, doHardCuts{false}, doEventCuts{false};
        bool doRotate{false};
        double Calculate(const std::vector<double>&, const double&, const int&);
        bool MakeCuts(Event&);
        bool MakeEventCuts(Event&);
        void Rotate(Event&);

        std::shared_ptr<Beam> beam;
        std::shared_ptr<Nucleus> nucleus;
        std::shared_ptr<Cascade> cascade;
        std::shared_ptr<HardScattering> scattering;
        Cuts hard_cuts{};
        Cuts event_cuts{};
        Vegas integrator;
        YAML::Node config;
        // std::vector<std::unique_ptr<HardScattering>> scatterings;
        // std::vector<Histogram> xsecsVsE;
        // std::vector<Vegas> integrators;
        
        std::shared_ptr<EventWriter> writer;

        Histogram hist;
};

}

#endif
