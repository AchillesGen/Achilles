#ifndef HEPMC3_EVENT_WRITER_HH
#define HEPMC3_EVENT_WRITER_HH

#include "Achilles/EventWriter.hh"
#include "Achilles/Statistics.hh"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wconversion"
#if defined(__clang__)
#pragma GCC diagnostic ignored "-Wdouble-promotion"
#elif defined(__GNUC__) || defined(__GNUG__)
#endif
#include "HepMC3/WriterAscii.h"
#pragma GCC diagnostic pop

namespace achilles {

class HepMC3Writer : public EventWriter {
  public:
    HepMC3Writer(const std::string& filename, bool zipped = true):
      file{HepMC3::WriterAscii(*EventWriter::InitializeStream(filename, zipped))}
      {}
    ~HepMC3Writer() override = default;

    void WriteHeader(const std::string &, const std::vector<ProcessGroup> &) override;
    void Write(const Event &) override;

  private:
    HepMC3::WriterAscii file;
    achilles::StatsData results;
};

} // namespace achilles

#endif
