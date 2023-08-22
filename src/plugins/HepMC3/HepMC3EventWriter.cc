#include "plugins/HepMC3/HepMC3EventWriter.hh"
#ifdef GZIP
#include "gzstream/gzstream.h"
#endif
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"

#include "Achilles/Event.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/Particle.hh"
#include "Achilles/Version.hh"
#include "spdlog/spdlog.h"

using achilles::HepMC3Writer;
using namespace HepMC3;

std::shared_ptr<std::ostream> HepMC3Writer::InitializeStream(const std::string &filename,
                                                             bool zipped) {
    std::shared_ptr<std::ostream> output = nullptr;
#ifdef GZIP
    if(zipped) {
        std::string zipname = filename;
        if(filename.substr(filename.size() - 3) != ".gz") zipname += std::string(".gz");
        output = std::make_shared<ogzstream>(zipname.c_str());
    } else
#endif
        output = std::make_shared<std::ofstream>(filename);

    return output;
}

void HepMC3Writer::WriteHeader(const std::string &filename, const std::vector<ProcessGroup> &) {
    // Setup generator information
    spdlog::trace("Writing Header");
    auto run = std::make_shared<GenRunInfo>();
    struct GenRunInfo::ToolInfo generator = {std::string("Achilles"), std::string(ACHILLES_VERSION),
                                             std::string("Neutrino event generator")};
    run->tools().push_back(generator);
    struct GenRunInfo::ToolInfo config = {std::string(filename), "1.0", std::string("Run card")};
    run->tools().push_back(config);
    run->set_weight_names({"Default"});

    file.set_run_info(run);
    file.write_run_info();
    spdlog::trace("Finished writing Header");
}

int ToHepMC3(achilles::ParticleStatus status) {
    switch(status) {
    case achilles::ParticleStatus::internal_test:
    case achilles::ParticleStatus::external_test:
    case achilles::ParticleStatus::propagating:
    case achilles::ParticleStatus::background:
    case achilles::ParticleStatus::captured:
        return 11;
    case achilles::ParticleStatus::initial_state:
        return 3;
    case achilles::ParticleStatus::final_state:
    case achilles::ParticleStatus::escaped:
        return 1;
    case achilles::ParticleStatus::decayed:
        return 2;
    case achilles::ParticleStatus::beam:
    case achilles::ParticleStatus::target:
        return 4;
    }
    return -1;
}

GenParticlePtr ToHepMC3(const achilles::Particle &particle) {
    HepMC3::FourVector mom{particle.Px(), particle.Py(), particle.Pz(), particle.E()};
    return std::make_shared<GenParticle>(mom, static_cast<int>(particle.ID()),
                                         ToHepMC3(particle.Status()));
}

struct HepMC3Visitor : achilles::HistoryVisitor {
    static constexpr double to_mm = 1e-12;
    GenEvent evt;
    struct compare {
        bool operator()(const achilles::Particle &a, const achilles::Particle &b) const {
            if(achilles::compare_momentum(a, 1e-10)(b)) { return a.Status() < b.Status(); }
            return a.Momentum().Vec3().Magnitude() < b.Momentum().Vec3().Magnitude();
        }
    };
    std::map<achilles::Particle, GenParticlePtr, HepMC3Visitor::compare> converted;
    std::vector<GenParticlePtr> beamparticles;
    HepMC3Visitor() : evt(Units::MEV, Units::MM), beamparticles(2) {}
    void visit(achilles::EventHistoryNode *node) {
        auto position = node->Position();
        HepMC3::FourVector vertex_pos{position.X(), position.Y(), position.Z(), 0};
        vertex_pos *= to_mm;
        GenVertexPtr vertex = std::make_shared<GenVertex>(vertex_pos);
        vertex->set_status(static_cast<int>(node->Status()));
        for(const auto &part : node->ParticlesIn()) {
            GenParticlePtr particle;
            if(converted.count(part) > 0) {
                particle = converted[part];
            } else {
                particle = ToHepMC3(part);
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
                particle = ToHepMC3(part);
                converted[part] = particle;
            }
            vertex->add_particle_out(particle);
        }
        evt.add_vertex(vertex);
    }
};

void HepMC3Writer::Write(const achilles::Event &event) {
    spdlog::debug("Writing out event");
    constexpr double to_mm = 1e-12;
    constexpr double nb_to_pb = 1000;

    // Update cumulative results, but skip writing if weight is zero
    results += event.Weight() * nb_to_pb;
    spdlog::trace("Event weight = {}", event.Weight());
    if(event.Weight() == 0) { return; }

    // Setup event units
    spdlog::trace("Setting up units");
    HepMC3Visitor visitor;
    visitor.evt.set_run_info(file.run_info());
    visitor.evt.set_event_number(results.Calls());

    // Interaction type
    // TODO: Add interaction type to the event, and have ids for different modes
    visitor.evt.add_attribute("InteractionType", std::make_shared<IntAttribute>(1));

    // Cross Section
    spdlog::trace("Writing out cross-section");
    auto cross_section = std::make_shared<GenCrossSection>();
    cross_section->set_cross_section(results.Mean(), results.Error(), results.FiniteCalls(),
                                     results.Calls());
    visitor.evt.add_attribute("GenCrossSection", cross_section);
    visitor.evt.add_attribute("Flux", std::make_shared<DoubleAttribute>(event.Flux()));
    visitor.evt.weight("Default") = event.Weight() * nb_to_pb;

    // TODO: once we have a detector to simulate interaction location
    // Event position
    // FourVector position{event.Position()};
    // evt.shift_position_to(position);

    // Walk the history and add to file
    event.History().WalkHistory(visitor);
    // visitor.evt.add_tree(visitor.beamparticles);
    file.write_event(visitor.evt);
}
