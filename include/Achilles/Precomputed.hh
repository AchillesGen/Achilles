#pragma once

#include "Achilles/Cascade.hh"
#include "Achilles/Event.hh"
#include "Achilles/EventWriter.hh"
#include "Achilles/Nucleus.hh"
#include "HepMC3/GenParticle_fwd.h"

#include <fstream>
#include <memory>
#include <optional>
#include <string>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wuseless-cast"
#include "NuHepMC/Reader.hxx"
#include "NuHepMC/ReaderUtils.hxx"
#include "NuHepMC/Types.hxx"
#pragma GCC diagnostic pop

namespace achilles {
class Event;
}

namespace achilles::Precomputed {

// Event Parsers
class EventReader {
  public:
    EventReader() = default;
    virtual ~EventReader() = default;

    virtual std::optional<Event> Next() = 0;
};

template <typename Derived>
using RegistrableEventReader =
    Registrable<EventReader, Derived, const std::string &, const std::shared_ptr<Nucleus> &>;
using EventReaderFactory =
    Factory<EventReader, const std::string &, const std::shared_ptr<Nucleus> &>;

class RawEventReader : public EventReader, RegistrableEventReader<RawEventReader> {
  public:
    RawEventReader(const std::string &, std::shared_ptr<Nucleus>);
    std::optional<Event> Next();

    // Required factory methods
    static std::unique_ptr<EventReader> Construct(const std::string &,
                                                  const std::shared_ptr<Nucleus> &);
    static std::string Name() { return "Raw"; }

  private:
    Event ParseLine(const std::string &);
    std::ifstream event_stream;
    std::shared_ptr<Nucleus> m_nucleus;

    double to_mev = 1;
};

class NuHepMCReader : public EventReader, RegistrableEventReader<NuHepMCReader> {
  public:
    NuHepMCReader(const std::string &, std::shared_ptr<Nucleus>);
    std::optional<Event> Next();

    // Required factory methods
    static std::unique_ptr<EventReader> Construct(const std::string &,
                                                  const std::shared_ptr<Nucleus> &);
    static std::string Name() { return "NuHepMC"; }

  private:
    Particle FromHepMC3(const HepMC3::ConstGenParticlePtr &) const;

    std::unique_ptr<NuHepMC::Reader> m_reader;
    std::shared_ptr<Nucleus> m_nucleus;
    NuHepMC::StatusCodeDescriptors m_vertex_status;
    NuHepMC::StatusCodeDescriptors m_particle_status;
};

void ConvertEvent(Event &, Particles &);

template <typename Reader, typename Pipeline> void RunPipeline(Reader &reader, Pipeline &pipeline) {
    while(auto evt = reader->Next()) { pipeline.Process(evt.value()); }
}

template <typename... Actions> class EventPipeline {
  public:
    EventPipeline(Actions... actions) : actions_(std::move(actions)...) {}

    void Process(Event &ev) {
        std::apply([&](auto &...action) { (action(ev), ...); }, actions_);
    }

  private:
    std::tuple<Actions...> actions_;
};

class RunCascade {
  public:
    RunCascade(const std::string &config);
    void RunAll();
    void Run(Event &event);

  private:
    Cascade cascade;
    std::shared_ptr<Nucleus> nucleus;
    std::unique_ptr<EventWriter> writer;
    std::unique_ptr<EventReader> reader;
};

} // namespace achilles::Precomputed
