#include "plugins/NuHepMC/NuHepMCWriter.hh"
#include "gzstream/gzstream.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/GenParticle.h"

#include "Achilles/Event.hh"
#include "Achilles/Version.hh"
#include "Achilles/Particle.hh"
#include "Achilles/Nucleus.hh"
#include "spdlog/spdlog.h"

using achilles::NuHepMCWriter;
using namespace HepMC3;

std::shared_ptr<std::ostream> NuHepMCWriter::InitializeStream(const std::string &filename, bool zipped) {
    std::shared_ptr<std::ostream> output = nullptr;
    if(zipped) {
        std::string zipname = filename;
        if(filename.substr(filename.size() - 3) != ".gz")
            zipname += std::string(".gz");
        output = std::make_shared<ogzstream>(zipname.c_str());
    } else {
        output = std::make_shared<std::ofstream>(filename);
    }

    return output;
}

void NuHepMCWriter::WriteHeader(const std::string &filename) {
    // Setup generator information
    spdlog::trace("Writing Header");
    auto run = std::make_shared<HepMC3::GenRunInfo>(); 
    run->add_attribute("NuHepMC.Version.Major",
                       std::make_shared<HepMC3::IntAttribute>(version[0]));
    run->add_attribute("NuHepMC.Version.Minor",
                       std::make_shared<HepMC3::IntAttribute>(version[1]));
    run->add_attribute("NuHepMC.Version.Patch",
                       std::make_shared<HepMC3::IntAttribute>(version[2]));
                                                  
    struct HepMC3::GenRunInfo::ToolInfo generator={std::string("Achilles"),
                                                   std::string(ACHILLES_VERSION),
                                                   std::string("Neutrino event generator")};
    run->tools().push_back(generator);
    run->add_attribute("Achilles.RunCard",
                       std::make_shared<HepMC3::StringAttribute>(filename));
    run->set_weight_names({"CV"});

    // Add all possible processes
    // TODO: Look up available processes and add bibtex info
    std::vector<int> proc_ids{101};
    run->add_attribute("NuHepMC.ProcessIDs",
                       std::make_shared<HepMC3::VectorIntAttribute>(proc_ids));
    run->add_attribute("NuHepMC.ProcessInfo[101].Name",
                       std::make_shared<HepMC3::StringAttribute>("SpectralQE"));
    run->add_attribute("NuHepMC.ProcessInfo[101].Description",
                       std::make_shared<HepMC3::StringAttribute>("Spectral function Quasielastic"));
    run->add_attribute("NuHepMC.ProcessInfo[101].Bibtex",
                       std::make_shared<HepMC3::StringAttribute>("Rocco:xxxx"));

    // List all possible vertex status codes
    // TODO: Make this a conversion from enum of the EventHistory class?
    std::vector<int> vertex_ids{1,2,3,4,5};
    run->add_attribute("NuHepMC.VertexStatusIDs",
                       std::make_shared<HepMC3::VectorIntAttribute>(vertex_ids));
    run->add_attribute("NuHepMC.VertexStatusInfo[1].Name",
                       std::make_shared<HepMC3::StringAttribute>("Primary"));
    run->add_attribute("NuHepMC.VertexStatusInfo[1].Description",
                       std::make_shared<HepMC3::StringAttribute>("The main hard interaction"));
    run->add_attribute("NuHepMC.VertexStatusInfo[2].Name",
                       std::make_shared<HepMC3::StringAttribute>("Beam"));
    run->add_attribute("NuHepMC.VertexStatusInfo[2].Description",
                       std::make_shared<HepMC3::StringAttribute>("The vertex defining a beam particle from a flux"));
    run->add_attribute("NuHepMC.VertexStatusInfo[3].Name",
                       std::make_shared<HepMC3::StringAttribute>("Nucleus"));
    run->add_attribute("NuHepMC.VertexStatusInfo[3].Description",
                       std::make_shared<HepMC3::StringAttribute>("The vertex defining a nucleon in a nucleus"));
    run->add_attribute("NuHepMC.VertexStatusInfo[4].Name",
                       std::make_shared<HepMC3::StringAttribute>("Decay"));
    run->add_attribute("NuHepMC.VertexStatusInfo[4].Description",
                       std::make_shared<HepMC3::StringAttribute>("A vertex used in a decay chain"));
    run->add_attribute("NuHepMC.VertexStatusInfo[5].Name",
                       std::make_shared<HepMC3::StringAttribute>("Shower"));
    run->add_attribute("NuHepMC.VertexStatusInfo[5].Description",
                       std::make_shared<HepMC3::StringAttribute>("Vertex used in a shower"));

    // List all possible particle status codes
    // TODO: Make this a conversion from ParticleStatus enum
    std::vector<int> particle_status{1, 2, 3, 4, 11};
    run->add_attribute("NuHepMC.ParticleStatusIDs",
                       std::make_shared<HepMC3::VectorIntAttribute>(particle_status));
    run->add_attribute("NuHepMC.ParticleStatusInfo[1].Name",
                       std::make_shared<HepMC3::StringAttribute>("Undecayed physical particle"));
    run->add_attribute("NuHepMC.ParticleStatusInfo[1].Description",
                       std::make_shared<HepMC3::StringAttribute>("Final state \"stable\" particles"));
    run->add_attribute("NuHepMC.ParticleStatusInfo[2].Name",
                       std::make_shared<HepMC3::StringAttribute>("Decayed physical particle"));
    run->add_attribute("NuHepMC.ParticleStatusInfo[2].Description",
                       std::make_shared<HepMC3::StringAttribute>("Particle that decayed during the generation"));
    run->add_attribute("NuHepMC.ParticleStatusInfo[3].Name",
                       std::make_shared<HepMC3::StringAttribute>("Documentation line"));
    run->add_attribute("NuHepMC.ParticleStatusInfo[3].Description",
                       std::make_shared<HepMC3::StringAttribute>("Internal particle history"));
    run->add_attribute("NuHepMC.ParticleStatusInfo[4].Name",
                       std::make_shared<HepMC3::StringAttribute>("Incoming beam Particle"));
    run->add_attribute("NuHepMC.ParticleStatusInfo[4].Description",
                       std::make_shared<HepMC3::StringAttribute>("Incoming beam particle"));
    run->add_attribute("NuHepMC.ParticleStatusInfo[11].Name",
                       std::make_shared<HepMC3::StringAttribute>("Target particle"));
    run->add_attribute("NuHepMC.ParticleStatusInfo[11].Description",
                       std::make_shared<HepMC3::StringAttribute>("Target particle"));

    // Signal conventions
    // TODO: Make flags to turn on / off different conventions
    std::vector<std::string> conventions{"G.C.2", "G.C.5", "E.C.1", "E.C.4", "E.C.5"};
    run->add_attribute("NuHepMC.Conventions",
                       std::make_shared<HepMC3::VectorStringAttribute>(conventions));

    // Write out the number of requested events
    // TODO: Read this from run card
    long nevents = 10;
    run->add_attribute("NuHepMC.Exposure.NEvents",
                       std::make_shared<HepMC3::LongAttribute>(nevents));

    file.set_run_info(run);
    file.write_run_info();
    spdlog::trace("Finished writing Header");
}

int ToNuHepMC(achilles::ParticleStatus status) {
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

GenParticlePtr ToNuHepMC(const achilles::Particle &particle) {
    HepMC3::FourVector mom{particle.Px(), particle.Py(), particle.Pz(), particle.E()};
    return std::make_shared<GenParticle>(mom, static_cast<int>(particle.ID()), ToNuHepMC(particle.Status()));
}

struct NuHepMCVisitor : achilles::HistoryVisitor {
    static constexpr double to_mm = 1e-12;
    GenEvent evt;
    struct compare {
        bool operator()(const achilles::Particle &a, const achilles::Particle &b) const {
            if(achilles::compare_momentum(a, 1e-10)(b)) {
                return a.Status() < b.Status();
            }
            return a.Momentum().Vec3().Magnitude() < b.Momentum().Vec3().Magnitude();
        }
    };
    std::map<achilles::Particle, GenParticlePtr, NuHepMCVisitor::compare> converted;
    std::vector<GenParticlePtr> beamparticles;
    NuHepMCVisitor() : evt(Units::MEV, Units::MM), beamparticles(2) {}
    void visit(achilles::EventHistoryNode *node) {
        auto position = node -> Position();
        HepMC3::FourVector vertex_pos{position.X(), position.Y(), position.Z(), 0};
        vertex_pos *= to_mm;
        GenVertexPtr vertex = std::make_shared<GenVertex>(vertex_pos);
        vertex->set_status(static_cast<int>(node->Status()));
        for(const auto &part : node -> ParticlesIn()) {
            GenParticlePtr particle;
            if(converted.count(part) > 0) {
                particle = converted[part];
            } else {
                particle = ToNuHepMC(part);
                converted[part] = particle;
            }
            if(node -> Status() == achilles::EventHistory::StatusCode::beam) {
                beamparticles[0] = particle;
            } else if (node -> Status() == achilles::EventHistory::StatusCode::target) {
                beamparticles[1] = particle;
            }
            vertex -> add_particle_in(particle);
        }
        for(const auto &part : node -> ParticlesOut()) {
            GenParticlePtr particle;
            if(converted.count(part) > 0) {
                particle = converted[part];
            } else {
                particle = ToNuHepMC(part);
                converted[part] = particle;
            }
            vertex -> add_particle_out(particle);
        }
        evt.add_vertex(vertex);
    }
};

void NuHepMCWriter::Write(const achilles::Event &event) {
    spdlog::debug("Writing out event");
    constexpr double to_mm = 1e-12;
    constexpr double nb_to_pb = 1000;

    // Update cumulative results, but skip writing if weight is zero
    results += event.Weight()*nb_to_pb;
    spdlog::trace("Event weight = {}", event.Weight());
    if(event.Weight() == 0) {
        return;
    }

    // Setup event units
    spdlog::trace("Setting up units");
    NuHepMCVisitor visitor;
    visitor.evt.set_run_info(file.run_info());
    visitor.evt.set_event_number(results.Calls());

    // Interaction type
    // TODO: Add interaction type to the event, and have ids for different modes
    visitor.evt.add_attribute("ProcID", std::make_shared<IntAttribute>(101));

    // Cross Section
    spdlog::trace("Writing out cross-section");
    auto cross_section = std::make_shared<GenCrossSection>();
    cross_section->set_cross_section(results.Mean(), results.Error(), results.FiniteCalls(), results.Calls());
    visitor.evt.add_attribute("GenCrossSection", cross_section);
    visitor.evt.add_attribute("Flux", std::make_shared<DoubleAttribute>(event.Flux()));
    visitor.evt.weight("CV") = event.Weight()*nb_to_pb;

    // TODO: once we have a detector to simulate interaction location
    // Event position
    // FourVector position{event.Position()};
    HepMC3::FourVector position{0, 0, 0, 0};
    visitor.evt.shift_position_to(position);
    std::vector<double> labpos = {0, 0, 0, 0};
    visitor.evt.add_attribute("LabPos", std::make_shared<VectorDoubleAttribute>(labpos));

    // Walk the history and add to file
    event.History().WalkHistory(visitor);
    // visitor.evt.add_tree(visitor.beamparticles);
    file.write_event(visitor.evt);
}
