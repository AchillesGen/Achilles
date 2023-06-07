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
    std::vector<int> vertex_ids{1,2};
    run->add_attribute("NuHepMC.VertexStatusIDs",
                       std::make_shared<HepMC3::VectorIntAttribute>(vertex_ids));
    run->add_attribute("NuHepMC.VertexStatusInfo[1].Name",
                       std::make_shared<HepMC3::StringAttribute>("Primary"));
    run->add_attribute("NuHepMC.VertexStatusInfo[1].Description",
                       std::make_shared<HepMC3::StringAttribute>("The main hard interaction"));
    run->add_attribute("NuHepMC.VertexStatusInfo[2].Name",
                       std::make_shared<HepMC3::StringAttribute>("Nucleus"));
    run->add_attribute("NuHepMC.VertexStatusInfo[2].Description",
                       std::make_shared<HepMC3::StringAttribute>("The vertex defining a nucleon in a nucleus"));

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

void NuHepMCWriter::Write(const achilles::Event &event) {
    using namespace HepMC3;
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
    evt.add_attribute("ProcID", std::make_shared<IntAttribute>(101));

    // Cross Section
    spdlog::trace("Writing out cross-section");
    auto cross_section = std::make_shared<GenCrossSection>();
    cross_section->set_cross_section(results.Mean(), results.Error(), results.FiniteCalls(), results.Calls());
    evt.add_attribute("GenCrossSection", cross_section);
    evt.add_attribute("Flux", std::make_shared<DoubleAttribute>(event.Flux()));
    evt.weight("CV") = event.Weight()*nb_to_pb;

    // TODO: once we have a detector to simulate interaction location
    // Event position
    // FourVector position{event.Position()};
    HepMC3::FourVector position{0, 0, 0, 0};
    evt.shift_position_to(position);
    std::vector<double> labpos = {0, 0, 0, 0};
    evt.add_attribute("LabPos", std::make_shared<VectorDoubleAttribute>(labpos));

    // Load in particle information
    const std::vector<achilles::Particle> hadrons = event.Hadrons();
    const std::vector<achilles::Particle> leptons = event.Leptons();
    // TODO: Clean this up
    const auto nuc_mass = achilles::ParticleInfo(event.CurrentNucleus() -> ID()).Mass();
    const HepMC3::FourVector initMass{0, 0, 0, nuc_mass};
    GenParticlePtr p1 = std::make_shared<GenParticle>(initMass, event.CurrentNucleus()->ID(), 4);
    HepMC3::FourVector hardVertexPos;
    GenParticlePtr nucleon;
    achilles::FourVector recoilMom{nuc_mass, 0, 0, 0};
    // TODO: Modify for MEC case
    size_t idx = 0;
    for(; idx < hadrons.size(); ++idx) {
        if(hadrons[idx].Status() == ParticleStatus::initial_state)
            break;
    }
    const auto initHadron = hadrons[idx];
    HepMC3::FourVector p2Mom{initHadron.Px(), initHadron.Py(), initHadron.Pz(), initHadron.E()};
    nucleon = std::make_shared<GenParticle>(p2Mom, int(initHadron.ID()), 3);

    // Add vertex for hadrons from nucleus
    const auto initPos = initHadron.Position();
    hardVertexPos = {initPos.X()*to_mm, initPos.Y()*to_mm, initPos.Z()*to_mm, 0};
    GenVertexPtr v1 = std::make_shared<GenVertex>(hardVertexPos);
    v1->set_status(2);
    v1->add_particle_in(p1);
    v1->add_particle_out(nucleon);
    evt.add_vertex(v1);

    // Add initial lepton
    // TODO: Get maximum neutrino energy from the beam
    //LP: Not sure about this 'Flux Choice Vertex', removed for now as it isn't in NuHepMC
    // const HepMC3::FourVector initBeam{0, 0, 10000, 10000};
    // GenParticlePtr p3In = std::make_shared<GenParticle>(initBeam, leptons[0].ID(), 4);
    const auto initLepton = leptons[0];
    const HepMC3::FourVector p3Mom{initLepton.Px(), initLepton.Py(), initLepton.Pz(), initLepton.E()};
    GenParticlePtr p3 = std::make_shared<GenParticle>(p3Mom, int(initLepton.ID()), 4);
    recoilMom += initLepton.Momentum();
    // Add vertex for neutrino flux
    // GenVertexPtr vLepIn = std::make_shared<GenVertex>(hardVertexPos);
    // vLepIn->add_particle_in(p3In);
    // vLepIn->add_particle_out(p3);
    // evt.add_vertex(vLepIn);

    // Add hard interaction vertex 
    GenVertexPtr v2 = std::make_shared<GenVertex>(hardVertexPos);
    v2->set_status(1);
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
    for(size_t i = 0; i < hadrons.size(); ++i) {
        if(i == idx) continue;
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
    if(event.Remnant().PID() != int(achilles::PID::dummyNucleus())) {
        const HepMC3::FourVector recoil{recoilMom.Px(), recoilMom.Py(), recoilMom.Pz(), recoilMom.E()};
        GenParticlePtr pRemnant = std::make_shared<GenParticle>(recoil, event.Remnant().PID(), 1);
        v2->add_particle_out(pRemnant);
    }
    evt.add_vertex(v2);

    file.write_event(evt);
}
