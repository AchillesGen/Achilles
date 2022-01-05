#include "plugins/HepMC3/HepMC3EventWriter.hh"
#include "gzstream/gzstream.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/GenParticle.h"

#include "nuchic/Event.hh"
#include "nuchic/Version.hh"
#include "nuchic/Particle.hh"
#include "nuchic/Nucleus.hh"
#include "spdlog/spdlog.h"

using nuchic::HepMC3Writer;
using namespace HepMC3;

std::shared_ptr<std::ostream> HepMC3Writer::InitializeStream(const std::string &filename, bool zipped) {
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

void HepMC3Writer::WriteHeader(const std::string &filename) {
    // Setup generator information
    spdlog::trace("Writing Header");
    auto run = std::make_shared<GenRunInfo>(); 
    struct GenRunInfo::ToolInfo generator={std::string("Nuchic"),
                                           std::string(NUCHIC_VERSION),
                                           std::string("Neutrino event generator")};
    run->tools().push_back(generator);
    struct GenRunInfo::ToolInfo config={std::string(filename),"1.0",std::string("Run card")};
    run->tools().push_back(config);
    run->set_weight_names({"Default"});

    file.set_run_info(run);
    file.write_run_info();
    spdlog::trace("Finished writing Header");
}

void HepMC3Writer::Write(const nuchic::Event &event) {
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
    GenEvent evt(Units::MEV, Units::MM);
    evt.set_run_info(file.run_info());
    evt.set_event_number(results.Calls());

    // Interaction type
    // TODO: Add interaction type to the event, and have ids for different modes
    evt.add_attribute("InteractionType", std::make_shared<IntAttribute>(1));

    // Cross Section
    spdlog::trace("Writing out cross-section");
    auto cross_section = std::make_shared<GenCrossSection>();
    cross_section->set_cross_section(results.Mean(), results.Error(), results.FiniteCalls(), results.Calls());
    evt.add_attribute("GenCrossSection", cross_section);
    evt.weight("Default") = event.Weight()*nb_to_pb;

    // TODO: once we have a detector to simulate interaction location
    // Event position
    // FourVector position{event.Position()};
    // evt.shift_position_to(position);

    // Load in particle information
    const std::vector<nuchic::Particle> hadrons = event.Hadrons();
    const std::vector<nuchic::Particle> leptons = event.Leptons();
    // TODO: Get nucleus mass from the nucleus object
    const HepMC3::FourVector initMass{0, 0, 0, 12000};
    GenParticlePtr p1 = std::make_shared<GenParticle>(initMass, event.CurrentNucleus()->ID(), 4);
    HepMC3::FourVector hardVertexPos;
    GenParticlePtr nucleon;
    nuchic::FourVector recoilMom{0, 0, 0, 12000};
    // TODO: Modify for MEC case
    const auto initHadron = hadrons[0];
    HepMC3::FourVector p2Mom{initHadron.Px(), initHadron.Py(), initHadron.Pz(), initHadron.E()};
    nucleon = std::make_shared<GenParticle>(p2Mom, int(initHadron.ID()), 3);

    // Add vertex for hadrons from nucleus
    const auto initPos = initHadron.Position();
    hardVertexPos = {initPos.X()*to_mm, initPos.Y()*to_mm, initPos.Z()*to_mm, 0};
    GenVertexPtr v1 = std::make_shared<GenVertex>(hardVertexPos);
    v1->add_particle_in(p1);
    v1->add_particle_out(nucleon);
    evt.add_vertex(v1);

    // Add initial lepton
    const auto initLepton = leptons[0];
    const HepMC3::FourVector p3Mom{initLepton.Px(), initLepton.Py(), initLepton.Pz(), initLepton.E()};
    GenParticlePtr p3 = std::make_shared<GenParticle>(p3Mom, int(initLepton.ID()), 4);
    recoilMom += initLepton.Momentum();

    // Add hard interaction vertex 
    GenVertexPtr v2 = std::make_shared<GenVertex>(hardVertexPos);
    v2->add_particle_in(nucleon);
    v2->add_particle_in(p3);

    // Add remaining leptons
    for(size_t i = 1; i < leptons.size(); ++i) {
        const auto currPart = leptons[i];
        const HepMC3::FourVector mom{currPart.Px(), currPart.Py(), currPart.Pz(), currPart.E()};
        GenParticlePtr p = std::make_shared<GenParticle>(mom, int(currPart.ID()), 1);
        v2->add_particle_out(p);
        recoilMom -= currPart.Momentum();
    }

    // Add remaining hard interaction hadrons
    for(size_t i = 1; i < hadrons.size(); ++i) {
        const auto currPart = hadrons[i];
        const HepMC3::FourVector mom{currPart.Px(), currPart.Py(), currPart.Pz(), currPart.E()};
        GenParticlePtr p = std::make_shared<GenParticle>(mom, int(currPart.ID()), 1);
        v2->add_particle_out(p);
        recoilMom -= currPart.Momentum();
    }

    // TODO: Add in cascade information if available

    // Add in remnant Nucleus
    // TODO: Get recoil momentum for the nucleus
    // TODO: Move remnant to last cascade vertex???
    const HepMC3::FourVector recoil{recoilMom.Px(), recoilMom.Py(), recoilMom.Pz(), recoilMom.E()};
    GenParticlePtr pRemnant = std::make_shared<GenParticle>(recoil, event.Remnant().PID(), 1);
    v2->add_particle_out(pRemnant);
    evt.add_vertex(v2);

    file.write_event(evt);
}
