#ifndef NUHEPMC_WRITER
#define NUHEPMC_WRITER

#include "Achilles/EventWriter.hh"
#include "Achilles/Statistics.hh"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wconversion"
#include "HepMC3/WriterAscii.h"
#pragma GCC diagnostic pop

namespace achilles {

class NuHepMCWriter : public EventWriter {
    public:
        NuHepMCWriter(const std::string &filename, bool zipped=true) 
            : file{InitializeStream(filename, zipped)} {}
        ~NuHepMCWriter() override = default;

        void WriteHeader(const std::string&) override;
        void Write(const Event&) override;

    private:
        static std::shared_ptr<std::ostream> InitializeStream(const std::string&, bool);
        HepMC3::WriterAscii file;
        achilles::StatsData results;
        static constexpr std::array<int, 3> version{0, 1, 0};
};

}

#endif
