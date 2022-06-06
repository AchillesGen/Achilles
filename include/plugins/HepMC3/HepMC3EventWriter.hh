#ifndef HEPMC3_EVENT_WRITER_HH
#define HEPMC3_EVENT_WRITER_HH

#include "Achilles/EventWriter.hh"
#include "Achilles/Statistics.hh"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wconversion"
#include "HepMC3/WriterAscii.h"
#pragma GCC diagnostic pop

namespace achilles {

class HepMC3Writer : public EventWriter {
    public:
        HepMC3Writer(const std::string &filename, bool zipped=true) 
            : file{InitializeStream(filename, zipped)} {}
        ~HepMC3Writer() override = default;

        void WriteHeader(const std::string&) override;
        void Write(const Event&) override;

    private:
        static std::shared_ptr<std::ostream> InitializeStream(const std::string&, bool);
        HepMC3::WriterAscii file;
        achilles::StatsData results;
};

}

#endif
