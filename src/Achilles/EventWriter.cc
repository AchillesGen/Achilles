#include "Achilles/EventWriter.hh"
#include "Achilles/Event.hh"
#include "Achilles/Particle.hh"
#include "Achilles/Version.hh"
#include "fmt/format.h"

std::ostream *achilles::EventWriter::InitializeStream(const std::string &filename, bool zip) {
#ifdef GZIP
    if(zip) {
        std::string zipname = filename;
        if(filename.substr(filename.size() - 3) != ".gz") zipname += std::string(".gz");
        return new ogzstream(zipname.c_str());
    }
#endif
    return new std::ofstream(filename);
}

void achilles::AchillesWriter::WriteHeader(const std::string &filename,
                                           const std::vector<ProcessGroup> &) {
    *m_out << fmt::format("Achilles Version: {}\n", ACHILLES_VERSION);
    *m_out << fmt::format("{0:-^40}\n\n", "");

    std::ifstream input(filename);
    std::string line;
    while(std::getline(input, line)) { *m_out << fmt::format("{}\n", line); }
    *m_out << fmt::format("{0:-^40}\n\n", "");
}

void achilles::AchillesWriter::Write(const Event &event) {
    *m_out << fmt::format("Event: {}\n", ++nEvents);
    *m_out << fmt::format("  Particles:\n");
    for(const auto &part : event.Particles()) { *m_out << fmt::format("  - {}\n", part); }
    *m_out << fmt::format("  - {}\n", event.Remnant());
    *m_out << fmt::format("  Weight: {}\n", event.Weight());
}
