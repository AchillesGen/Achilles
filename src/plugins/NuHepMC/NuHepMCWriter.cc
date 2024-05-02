#include "plugins/NuHepMC/NuHepMCWriter.hh"

#include "Achilles/Event.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/Particle.hh"
#include "Achilles/Process.hh"
#include "Achilles/Version.hh"

#include "gzstream/gzstream.h"

#include "NuHepMC/WriterUtils.hxx"
#include "NuHepMC/make_writer.hxx"

#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"

#include "spdlog/spdlog.h"

using achilles::NuHepMCWriter;
using namespace HepMC3;

namespace NuHepMC {
namespace VertexStatus {
static constexpr int Propagation = 22;
static constexpr int Cascade = 23;
static constexpr int Beam = 24;
static constexpr int Decay = 25;
static constexpr int Shower = 26;
} // namespace VertexStatus

namespace ParticleStatus {
static constexpr int InternalTest = 22;
static constexpr int ExternalTest = 23;
static constexpr int Propagating = 24;
static constexpr int Background = 25;
static constexpr int Captured = 26;
static constexpr int UndecayedResidue = 27;
static constexpr int Spectator = 28;
static constexpr int Cascade = 29;
} // namespace ParticleStatus
} // namespace NuHepMC

NuHepMCWriter::NuHepMCWriter(const std::string &filename, bool zipped) : outfilename(filename) {}

void NuHepMCWriter::WriteHeader(const std::string &filename,
                                const std::vector<ProcessGroup> &groups) {
    // Setup generator information
    spdlog::trace("Writing Header");
    auto run = std::make_shared<HepMC3::GenRunInfo>();
    NuHepMC::GR2::WriteVersion(run);

    struct HepMC3::GenRunInfo::ToolInfo generator = {std::string("Achilles"),
                                                     std::string(ACHILLES_VERSION),
                                                     std::string("Neutrino event generator")};
    run->tools().push_back(generator);
    NuHepMC::add_attribute(run, "Achilles.RunCard", filename);

    NuHepMC::GR7::SetWeightNames(run, {
                                          "CV",
                                      });

    // Add all possible processes
    std::vector<int> proc_ids = achilles::AllProcessIDs(groups);
    run->add_attribute("NuHepMC.ProcessIDs",
                       std::make_shared<HepMC3::VectorIntAttribute>(proc_ids));
    auto metadata = achilles::AllProcessMetadata(groups);
    for(const auto data : metadata) {
        NuHepMC::add_attribute(run, fmt::format("NuHepMC.ProcessInfo[{}].Name", data.id),
                               data.name);
        NuHepMC::add_attribute(run, fmt::format("NuHepMC.ProcessInfo[{}].Description", data.id),
                               data.description);
        NuHepMC::add_attribute(run, fmt::format("NuHepMC.ProcessInfo[{}].InspireHEP", data.id),
                               data.inspireHEP);
    }

    // List all possible vertex status codes
    // TODO: Make this a conversion from enum of the EventHistory class?
    NuHepMC::GR5::WriteVertexStatusIDDefinitions(
        run, {
                 {NuHepMC::VertexStatus::Primary, {"Primary", "The main hard interaction"}},
                 {NuHepMC::VertexStatus::Beam,
                  {"Beam", "The vertex defining a beam particle from a flux"}},
                 {NuHepMC::VertexStatus::NucleonSeparation,
                  {"Nucleus", "The vertex defining a nucleon in a nucleus"}},
                 {NuHepMC::VertexStatus::Decay, {"Decay", "A vertex used in a decay chain"}},
                 {NuHepMC::VertexStatus::Shower, {"Shower", "Vertex used in a shower"}},
                 {NuHepMC::VertexStatus::Propagation, {"Propagation", "Propagation"}},
                 {NuHepMC::VertexStatus::Cascade, {"Cascade", "Cascade"}},
             });

    // List all possible particle status codes
    // TODO: Make this a conversion from ParticleStatus enum
    NuHepMC::GR6::WriteParticleStatusIDDefinitions(
        run, {
                 {NuHepMC::ParticleStatus::UndecayedPhysical,
                  {"Undecayed physical particle", "Final state \"stable\" particles"}},
                 {NuHepMC::ParticleStatus::DecayedPhysical,
                  {"Decayed physical particle", "Particle that decayed during the generation"}},
                 {NuHepMC::ParticleStatus::DocumentationLine,
                  {"Documentation line", "Internal particle history"}},
                 {NuHepMC::ParticleStatus::IncomingBeam,
                  {"Incoming beam Particle", "Incoming beam Particle"}},
                 {NuHepMC::ParticleStatus::Target, {"Target particle", "Target particle"}},
                 {NuHepMC::ParticleStatus::InternalTest, {"InternalTest", "InternalTest"}},
                 {NuHepMC::ParticleStatus::ExternalTest, {"ExternalTest", "ExternalTest"}},
                 {NuHepMC::ParticleStatus::Propagating, {"Propagating", "Propagating"}},
                 {NuHepMC::ParticleStatus::Background, {"Background", "Background"}},
                 {NuHepMC::ParticleStatus::Captured, {"Captured", "Captured"}},
                 {NuHepMC::ParticleStatus::Spectator, {"Spectator", "Spectator"}},
                 {NuHepMC::ParticleStatus::Cascade, {"Cascade", "Cascade"}},
             });

    // Signal conventions
    // TODO: Make flags to turn on / off different conventions
    NuHepMC::GC1::SetConventions(run,
                                 {"G.C.2", "G.C.4", "G.C.6", "E.C.1", "E.C.4", "E.C.5", "V.C.1"});

    NuHepMC::GC4::SetCrossSectionUnits(run, "pb", "PerTargetAtom");

    // Write out the number of requested events
    // TODO: Read this from run card
    long nevents = 10;
    NuHepMC::GC2::SetExposureNEvents(run, nevents);

    file = std::shared_ptr<HepMC3::Writer>(NuHepMC::Writer::make_writer(outfilename, run));
    spdlog::trace("Finished writing Header");
}

int ToNuHepMC(achilles::ParticleStatus status) {
    switch(status) {
    case achilles::ParticleStatus::internal_test:
        return NuHepMC::ParticleStatus::InternalTest;
    case achilles::ParticleStatus::external_test:
        return NuHepMC::ParticleStatus::ExternalTest;
    case achilles::ParticleStatus::propagating:
        return NuHepMC::ParticleStatus::Propagating;
    case achilles::ParticleStatus::background:
        return NuHepMC::ParticleStatus::Background;
    case achilles::ParticleStatus::captured:
        return NuHepMC::ParticleStatus::Captured;
    case achilles::ParticleStatus::initial_state:
        return NuHepMC::ParticleStatus::DecayedPhysical;
    case achilles::ParticleStatus::final_state:
    case achilles::ParticleStatus::escaped:
        return NuHepMC::ParticleStatus::UndecayedPhysical;
    case achilles::ParticleStatus::decayed:
        return NuHepMC::ParticleStatus::DecayedPhysical;
    case achilles::ParticleStatus::beam:
        return NuHepMC::ParticleStatus::IncomingBeam;
    case achilles::ParticleStatus::target:
        return NuHepMC::ParticleStatus::Target;
    case achilles::ParticleStatus::spectator:
        return NuHepMC::ParticleStatus::Spectator;
    case achilles::ParticleStatus::interacted:
        return NuHepMC::ParticleStatus::Cascade;
    }
    return -1;
}

int ToNuHepMC(achilles::EventHistoryNode::StatusCode status) {
    switch(status) {
    case achilles::EventHistoryNode::StatusCode::propagation:
        return NuHepMC::VertexStatus::Propagation;
    case achilles::EventHistoryNode::StatusCode::cascade:
        return NuHepMC::VertexStatus::Cascade;
    case achilles::EventHistoryNode::StatusCode::primary:
        return NuHepMC::VertexStatus::Primary;
    case achilles::EventHistoryNode::StatusCode::beam:
        return NuHepMC::VertexStatus::Beam;
    case achilles::EventHistoryNode::StatusCode::target:
        return NuHepMC::VertexStatus::NucleonSeparation;
    case achilles::EventHistoryNode::StatusCode::decay:
        return NuHepMC::VertexStatus::Decay;
    case achilles::EventHistoryNode::StatusCode::shower:
        return NuHepMC::VertexStatus::Shower;
    }
    return -1;
}

GenParticlePtr ToNuHepMC(const achilles::Particle &particle) {
    HepMC3::FourVector mom{particle.Px(), particle.Py(), particle.Pz(), particle.E()};
    return std::make_shared<GenParticle>(mom, static_cast<int>(particle.ID()),
                                         ToNuHepMC(particle.Status()));
}

struct NuHepMCVisitor : achilles::HistoryVisitor {
    static constexpr double to_mm = 1e-12;
    GenEvent evt;
    struct compare {
        bool operator()(const achilles::Particle &a, const achilles::Particle &b) const {
            if(a.ID() != b.ID()) return a.ID() < b.ID();
            if(achilles::compare_momentum(a, 1e-10)(b)) { return a.Status() < b.Status(); }
            return a.Momentum().Vec3().Magnitude() < b.Momentum().Vec3().Magnitude();
        }
    };
    std::map<achilles::Particle, GenParticlePtr, NuHepMCVisitor::compare> converted;
    std::vector<GenParticlePtr> beamparticles;
    NuHepMCVisitor() : evt(Units::MEV, Units::MM), beamparticles(2) {}
    void visit(achilles::EventHistoryNode *node) {
        auto position = node->Position();
        HepMC3::FourVector vertex_pos{position.X(), position.Y(), position.Z(), 0};
        vertex_pos *= to_mm;
        GenVertexPtr vertex = std::make_shared<GenVertex>(vertex_pos);
        vertex->set_status(ToNuHepMC(node->Status()));
        for(const auto &part : node->ParticlesIn()) {
            GenParticlePtr particle;
            if(converted.count(part) > 0) {
                particle = converted[part];
            } else {
                particle = ToNuHepMC(part);
                converted[part] = particle;
            }
            if(node->Status() == achilles::EventHistory::StatusCode::beam) {
                beamparticles[0] = particle;
            } else if(node->Status() == achilles::EventHistory::StatusCode::target) {
                beamparticles[1] = particle;
            }
            vertex->add_particle_in(particle);
        }
        for(const auto &part : node->ParticlesOut()) {
            GenParticlePtr particle;
            if(converted.count(part) > 0) {
                particle = converted[part];
            } else {
                particle = ToNuHepMC(part);
                converted[part] = particle;
            }
            vertex->add_particle_out(particle);
        }
        evt.add_vertex(vertex);
    }
};

void NuHepMCWriter::Write(const achilles::Event &event) {
    spdlog::debug("Writing out event");
    constexpr double to_mm = 1e-12;
    constexpr double nb_to_pb = 1000;

    // Update cumulative results, but skip writing if weight is zero
    results += event.Weight() * nb_to_pb;
    spdlog::trace("Event weight = {}", event.Weight());
    if(event.Weight() == 0) { return; }

    // Setup event units
    spdlog::trace("Setting up units");
    NuHepMCVisitor visitor;
    visitor.evt.set_run_info(file->run_info());
    visitor.evt.set_event_number(results.Calls());

    // Interaction type
    NuHepMC::ER3::SetProcessID(visitor.evt, event.ProcessId());

    // Cross Section
    spdlog::trace("Writing out cross-section");
    auto cross_section = std::make_shared<GenCrossSection>();
    cross_section->set_cross_section(results.Mean(), results.Error(), results.FiniteCalls(),
                                     results.Calls());
    visitor.evt.set_cross_section(cross_section);

    NuHepMC::add_attribute(visitor.evt, "Flux", event.Flux());
    visitor.evt.weight("CV") = event.Weight() * nb_to_pb;

    // TODO: once we have a detector to simulate interaction location
    // Event position
    // FourVector position{event.Position()};
    HepMC3::FourVector position{0, 0, 0, 0};
    visitor.evt.shift_position_to(position);
    NuHepMC::ER5::SetLabPosition(visitor.evt, {0, 0, 0, 0});

    // Walk the history and add to file
    event.History().WalkHistory(visitor);
    // visitor.evt.add_tree(visitor.beamparticles);
    file->write_event(visitor.evt);
}
