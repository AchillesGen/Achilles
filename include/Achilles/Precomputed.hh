#pragma once

#include "Achilles/Cascade.hh"
#include "Achilles/Event.hh"
#include "Achilles/EventWriter.hh"
#include "Achilles/Nucleus.hh"

#include <fstream>
#include <memory>
#include <optional>
#include <string>

namespace achilles {
class Event;
}

namespace achilles::Precomputed {

// Event Parsers
class RawEventReader {
  public:
    RawEventReader(const std::string &, std::shared_ptr<Nucleus>);
    std::optional<Event> Next();

  private:
    Event ParseLine(const std::string &);
    std::ifstream event_stream;
    std::shared_ptr<Nucleus> m_nucleus;

    double to_mev = 1;
};

class NuHepMCReader {
  public:
    NuHepMCReader(const std::string &, std::shared_ptr<Nucleus>);
    std::optional<Event> Next();

  private:
    std::shared_ptr<Nucleus> m_nucleus;
};

void ConvertEvent(Event &, Particles &);

template <typename Reader, typename Pipeline> void RunPipeline(Reader &reader, Pipeline &pipeline) {
    while(auto evt = reader.Next()) { pipeline.Process(evt.value()); }
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
    std::string event_filename;
};

} // namespace achilles::Precomputed
