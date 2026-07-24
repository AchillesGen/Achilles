#ifndef EVENTGEN_HH
#define EVENTGEN_HH

#include "Achilles/CombinedCuts.hh"
#include "Achilles/Event.hh"
#include "Achilles/Histogram.hh"
#include "Achilles/MultiChannel.hh"
#include "Achilles/ParticleInfo.hh"
#include "Achilles/Process.hh"
#include "Achilles/QuasielasticTestMapper.hh"
#include "Achilles/Settings.hh"
#include "Achilles/Unweighter.hh"
#include "Achilles/Vegas.hh"

#include <memory>
#include <optional>
#include <utility>
#include <vector>

namespace YAML {
class Node;
}

namespace HepMC3 {
class GenRunInfo;
}

namespace achilles {

// Forward declare types
class Event;
class Beam;
class GeometryBeam;
class FourVector;
class Nucleus;
class Cascade;
class EventWriter;

class SherpaInterface;

class EventGen {
  public:
    EventGen(const std::string &, std::vector<std::string>);
    // Construct from an already-loaded (possibly programmatically updated)
    // YAML node instead of a file on disk.  Nested !include tags are still
    // resolved.  `run_card_name` is only recorded as the Achilles.RunCard
    // run-info attribute for provenance.
    EventGen(const YAML::Node &node, std::vector<std::string> shargs,
             const std::string &run_card_name = "<YAML node>");
    void Initialize();
    void GenerateEvents(bool);

    // Geometry-driver API (used by external geometry codes such as NuGeometry).
    Beam *GetBeam() { return beam.get(); }
    // Populate an externally owned NuHepMC run_info (e.g. NuGeometry's) with the
    // Achilles GR.* run-level metadata. Requires HepMC3 support.
    void WriteRunInfo(std::shared_ptr<HepMC3::GenRunInfo> run) const;
    // Inject one lab-frame ray of species nu on the given target; restricts the
    // next GenerateSingleEvent() to that (nu, target).
    void InjectRay(const FourVector &lab_ray, PID nu, PID target);
    bool GenerateSingleEvent();
    // The most recent event from GenerateSingleEvent(), in the lab frame (the
    // geometry rotation that aligned the injected ray with +z has been undone).
    // Valid only after GenerateSingleEvent() returned true.
    const Event &LastEvent() const { return m_last_event; }
    void Finalize();

    // Unweighting envelope (max-weight, in units of sigma) for the matching
    // processes of an injected (nu, target); the accept/reject ceiling NuGeom
    // uses in EnvelopeNoRetry mode. Valid after Initialize().
    double VertexEnvelope(PID nu, PID target) const;

    // Total cross section sigma(E) for an injected (nu, target), summed over all
    // matching process groups on that target. This is NuGeom's TotalXSecRetry
    // layer-1 accessor. Valid after Initialize() with an energy-range beam.
    double TotalXSec(double energy, PID nu, PID target) const;

    // Which layer-1 envelope the geometry driver uses.  The physics path and
    // the per-event weight (partial-unweight max(1, raw_w/w_max_p) ~ 1) are
    // identical in both modes; the modes differ only in what the driver uses
    // for its layer-1 acceptance (VertexEnvelope vs the sigma_tot spline) and
    // whether it retries a vertex until emit.  The sigma scale for the
    // driver's estimator comes from VertexEnvelope (EnvelopeNoRetry) or
    // LastProcessMaxWeight (TotalXSecRetry).
    enum class SamplingMode { EnvelopeNoRetry, TotalXSecRetry };
    void SetSamplingMode(SamplingMode mode) { m_sampling_mode = mode; }
    SamplingMode GetSamplingMode() const { return m_sampling_mode; }

  private:
    // Shared constructor body (both the file-path and YAML-node constructors
    // delegate here once `config` and `m_run_card` are initialized).
    void Setup(std::vector<std::string> shargs);

    bool runCascade{false}, outputEvents{false}, doHardCuts{false};
    bool runDecays{false};
    bool doRotate{false};
    bool MakeCuts(Event &);
    // bool MakeEventCuts(Event&);
    void Rotate(Event &);

    std::shared_ptr<Beam> beam;
    std::vector<std::shared_ptr<Nucleus>> nuclei;
    std::shared_ptr<Cascade> cascade;
    std::vector<ProcessGroup> process_groups;
    CutCollection hard_cuts{};
    // CutCollection event_cuts{};
    MultiChannel integrator;
    Integrand<FourVector> integrand;
    Settings config;
    std::vector<double> m_group_weights{};
    double m_max_weight{};

    std::shared_ptr<EventWriter> writer;
    std::unique_ptr<Unweighter> unweighter;
    SherpaInterface *p_sherpa;

    // Geometry mode: non-owning pointer to the GeometryBeam (if the configured
    // beam is one), and the current injected (neutrino, target) for selection.
    GeometryBeam *m_geometry_beam{nullptr};
    std::optional<std::pair<PID, PID>> m_injection{};
    // Rate-limit the injected-energy-outside-EnergyRange warning.
    static constexpr size_t kMaxRangeWarnings = 10;
    size_t m_range_warnings{0};
    SamplingMode m_sampling_mode{SamplingMode::EnvelopeNoRetry};
    Event m_last_event{};
    std::string m_run_card{};
};

} // namespace achilles

#endif
