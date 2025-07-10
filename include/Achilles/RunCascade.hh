#pragma once

#include "Achilles/Cascade.hh"
#include "Achilles/EventWriter.hh"
#include "Achilles/Histogram.hh"
#include "Achilles/Settings.hh"
#include "yaml-cpp/yaml.h"

#include <string>

namespace achilles {

class Nucleus;

enum class CascadeMode {
    CrossSection,
    Transparency,
};

inline std::string ToString(CascadeMode mode) {
    switch(mode) {
    case CascadeMode::CrossSection:
        return "CrossSection";
    case CascadeMode::Transparency:
        return "Transparency";
    }
    return "Unknown";
}

} // namespace achilles

namespace YAML {

template <> struct convert<achilles::CascadeMode> {
    static bool decode(const YAML::Node &node, achilles::CascadeMode &mode) {
        std::string name = node.as<std::string>();
        if(name == "CrossSection")
            mode = achilles::CascadeMode::CrossSection;
        else if(name == "Transparency")
            mode = achilles::CascadeMode::Transparency;
        else
            return false;

        return true;
    }
};

} // namespace YAML

namespace achilles::CascadeTest {

void RunCascade(const std::string &runcard);

void InitCrossSection(Event &, PID, double, double, std::shared_ptr<Nucleus>);
void InitTransparency(Event &, PID, double, std::shared_ptr<Nucleus>, bool external = false);

class CascadeRunner {
  public:
    CascadeRunner(const std::string &);
    void run();

  private:
    void GenerateEvent(double, Histogram &, Histogram &);
    bool NeedsEvents() const { return generated_events < requested_events; }
    void Reset() { generated_events = 0; }

    size_t requested_events{}, generated_events{};

    std::map<std::string, double> m_params;
    std::string m_output_name;
    PID m_pid;

    CascadeMode m_mode;
    std::unique_ptr<Cascade> m_cascade;
    std::pair<double, double> m_mom_range;

    std::shared_ptr<Nucleus> m_nuc;
    std::unique_ptr<EventWriter> m_writer{nullptr};

    Settings setting;
};

} // namespace achilles::CascadeTest
