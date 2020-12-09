#ifndef EVENTGEN_HH
#define EVENTGEN_HH

#include "nuchic/Histogram.hh"
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

class EventGen {
    public:
        EventGen(const std::string&);
        void Initialize();
        void GenerateEvents();

    private:
        bool runCascade, outputEvents;
        double Calculate(const std::vector<double>&, const double&);

        std::shared_ptr<Beam> beam;
        std::shared_ptr<Nucleus> nucleus;
        std::shared_ptr<Cascade> cascade;
        std::shared_ptr<HardScattering> scattering;
        Vegas integrator;
        YAML::Node config;
        // std::vector<std::unique_ptr<HardScattering>> scatterings;
        // std::vector<Histogram> xsecsVsE;
        // std::vector<Vegas> integrators;
        
        std::shared_ptr<EventWriter> writer;
        std::shared_ptr<randutils::mt19937_rng> rng;

        Histogram hist;
};

}

#endif
