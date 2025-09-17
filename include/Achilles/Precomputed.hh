#pragma once

#include "Achilles/Cascade.hh"
#include "Achilles/Event.hh"
#include "Achilles/EventWriter.hh"
#include "Achilles/Nucleus.hh"

#include <fstream>
#include <memory>
#include <string>
#include <vector>

namespace achilles {
class Event;
}

namespace achilles::Precomputed {

Event ParseEvent(const std::string &line, std::shared_ptr<Nucleus> nuc);
template <typename Func>
void ProcessEventsFromFile(Func func, const std::string &event_file, std::shared_ptr<Nucleus> nuc) {
    std::ifstream ievents(event_file.c_str());
    std::string line;
    while(std::getline(ievents, line)) {
        auto event = ParseEvent(line, nuc);
        func(event);
    }
}

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
