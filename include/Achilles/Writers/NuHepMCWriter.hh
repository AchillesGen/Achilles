#ifndef NUHEPMC_WRITER
#define NUHEPMC_WRITER

#include "Achilles/Statistics.hh"
#include "Achilles/Writers/EventWriter.hh"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wconversion"
#include "HepMC3/Writer.h"
#pragma GCC diagnostic pop

#include <memory>

namespace HepMC3 {
class GenRunInfo;
class GenEvent;
} // namespace HepMC3

namespace achilles {

class ProcessGroup;
class Event;

// Merge a generated achilles::Event into an externally owned GenEvent (e.g. one
// NuGeometry pre-filled with the primary vertex, incoming neutrino, and target).
// The full Achilles history (primary outgoing + cascade) is appended onto the
// existing primary vertex, reusing NuGeom's incoming nu/target; Achilles values
// are converted into the event's units.
//
// `results` is the caller-owned running sigma estimator stamped as E.C.4; the
// caller must accumulate it *before* calling -- one entry per generator trial
// (zero for rejected trials, the sigma-units trial weight in pb for emitted
// ones) so that Mean() = sum(w)/N_trials stays an unbiased sigma estimate.
// The per-event weight is stamped as E.C.1 ("CV").
void FillHepMC3Event(const Event &event, HepMC3::GenEvent &evt, const StatsData &results);

class NuHepMCWriter : public EventWriter {
  public:
    NuHepMCWriter(const string &filename, bool) : outfilename{filename} {}
    ~NuHepMCWriter() override = default;

    void WriteHeader(const std::string &, const std::vector<ProcessGroup> &) override;
    // Populate a (possibly externally owned, e.g. NuGeometry's) run_info with the
    // Achilles GR.* run-level metadata. The file-writing WriteHeader override
    // delegates here. runcard is the path recorded as Achilles.RunCard.
    static void WriteHeader(std::shared_ptr<HepMC3::GenRunInfo> run, const std::string &runcard,
                            const std::vector<ProcessGroup> &groups);
    void Write(const Event &) override;

  private:
    std::shared_ptr<HepMC3::Writer> file;
    std::string outfilename;
    achilles::StatsData results;
    static constexpr std::array<int, 3> version{0, 1, 0};
};

} // namespace achilles

#endif
