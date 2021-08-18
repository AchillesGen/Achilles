#ifndef EVENTGEN_HH
#define EVENTGEN_HH

#include "nuchic/CombinedCuts.hh"
#include "nuchic/Histogram.hh"
#include "nuchic/ParticleInfo.hh"
#include "nuchic/QuasielasticTestMapper.hh"
#include "nuchic/Vegas.hh"
#include "nuchic/MultiChannel.hh"

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

class SherpaMEs;

class EventGen {
    public:
        EventGen(const std::string&,SherpaMEs *const);
        void Initialize();
        void GenerateEvents();
        size_t nevents;
        size_t total_events;
        size_t max_batch;

    private:
        bool runCascade, outputEvents, doHardCuts{false}, doEventCuts{false};
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
        SherpaMEs *sherpa;
        YAML::Node config;
        // std::vector<std::unique_ptr<HardScattering>> scatterings;
        // std::vector<Histogram> xsecsVsE;
        // std::vector<Vegas> integrators;

        std::shared_ptr<EventWriter> writer;

        Histogram hist, hist2, hist3, hist4, hist5, hist6;
};

}

#endif
