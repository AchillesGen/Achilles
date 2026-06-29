#include "Achilles/Precomputed.hh"
#include "Achilles/Cascade.hh"
#include "Achilles/EventHistory.hh"
#include "Achilles/EventWriter.hh"
#include "Achilles/Exception.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/Particle.hh"
#include "Achilles/Potential.hh"
#include "Achilles/Process.hh"
#include "Achilles/Random.hh"
#include "Achilles/Settings.hh"
#include "Achilles/Utilities.hh"
#include <functional>
#include <sstream>
#include <stdexcept>

#ifdef ACHILLES_ENABLE_HEPMC3
#include "HepMC3/Print.h"
#include "NuHepMC/EventUtils.hxx"
#include "Plugins/HepMC3/HepMC3EventWriter.hh"
#include "Plugins/NuHepMC/NuHepMCWriter.hh"
#endif

#include "spdlog/spdlog.h"
#include "yaml-cpp/yaml.h"

#include <chrono>

using namespace achilles::Precomputed;

bool iequals(const std::string &a, const std::string &b) {
    return a.size() == b.size() &&
           std::equal(a.begin(), a.end(), b.begin(), [](unsigned char c1, unsigned char c2) {
               return std::tolower(c1) == std::tolower(c2);
           });
}

std::unique_ptr<EventReader> RawEventReader::Construct(const std::string &filename,
                                                       const std::shared_ptr<Nucleus> &nuc) {
    return std::make_unique<RawEventReader>(filename, nuc);
}

RawEventReader::RawEventReader(const std::string &filename, std::shared_ptr<Nucleus> nuc)
    : event_stream(filename), m_nucleus(nuc) {
    if(!event_stream.is_open()) {
        auto err_msg = fmt::format("Precomputed: File {} does not exist!", filename);
        throw std::runtime_error(err_msg);
    }

    // Read the preamble
    // NOTE: Currently just the units exist in preamble
    std::string line;
    std::getline(event_stream, line);
    if(iequals("gev", line))
        to_mev = 1000;
    else if(!iequals("mev", line)) {
        throw std::runtime_error("Precomputed: Invalid unit found in preamble");
    }
}

std::optional<achilles::Event> RawEventReader::Next() {
    // Read the line for the event
    std::string line;
    if(!std::getline(event_stream, line)) return std::nullopt;

    return ParseLine(line);
}

achilles::Event RawEventReader::ParseLine(const std::string &line) {
    std::istringstream ss(line);

    // Read in event weight
    double wgt;
    ss >> wgt;

    // Read in all the particles
    std::vector<Particle> particles;
    while(ss) {
        long int pdg;
        double E, px, py, pz;
        ss >> pdg >> E >> px >> py >> pz;
        FourVector mom{E * to_mev, px * to_mev, py * to_mev, pz * to_mev};
        particles.emplace_back(pdg, mom);
    }

    // Setup event
    Event event(m_nucleus, {}, 1);
    event.Weight() = wgt;
    ConvertEvent(event, particles);

    return event;
}

void achilles::Precomputed::ConvertEvent(Event &event, Particles &particles) {
// Set the status of the particles and add to the event
    
	// MEC mode: the first two hadrons in the raw event are interpreted
	// as the two initial nucleons participating in the two-body current.
	std::size_t n_leptons_seen = 0;
	std::size_t n_initial_hadrons_seen = 0;
	std::size_t n_outgoing_hadrons_seen = 0;

	std::vector<ThreeVector> initial_positions;
	
	for(auto &particle : particles) {

	   if(particle.Info().IsLepton()) {
    	       if(n_leptons_seen == 0) {
        	   particle.Status() = ParticleStatus::initial_state;
               } else {
                  particle.Status() = ParticleStatus::final_state;
               }
               event.Leptons().push_back(particle);
               ++n_leptons_seen;

	    } else if(particle.Info().IsHadron()) {

		if(n_leptons_seen == 0) {
                   throw std::runtime_error(
                   "Precomputed: Found a hadron before the initial lepton");
                }    
		if(n_leptons_seen == 1) {
		    if(particle.ID() == PID::proton()) {
			auto protons = event.Protons(ParticleStatus::background);

			if(protons.empty()) {
			    throw std::runtime_error("Precomputed: No background protons available");
			}

			auto sampled_protons = Random::Instance().Sample(1, protons);
			auto sampled_proton = sampled_protons[0];

			// Use the momentum from the precomputed event,
			// but keep the position from the sampled Achilles nucleon.
			sampled_proton.get().Momentum() = particle.Momentum();
			sampled_proton.get().Status() = ParticleStatus::initial_state;

			initial_positions.push_back(sampled_proton.get().Position());

		    } else if(particle.ID() == PID::neutron()) {
			auto neutrons = event.Neutrons(ParticleStatus::background);

			if(neutrons.empty()) {
			    throw std::runtime_error("Precomputed: No background neutrons available");
			}

			auto sampled_neutrons = Random::Instance().Sample(1, neutrons);
			auto sampled_neutron = sampled_neutrons[0];

			sampled_neutron.get().Momentum() = particle.Momentum();
			sampled_neutron.get().Status() = ParticleStatus::initial_state;

			initial_positions.push_back(sampled_neutron.get().Position());

		    } else {
			throw std::runtime_error(
			    "Precomputed: Initial hadron must be a proton or neutron");
		    }

		    ++n_initial_hadrons_seen;

		} else {
		    
		    particle.Status() = ParticleStatus::propagating;

		    if(initial_positions.empty()) {
			throw std::runtime_error(
			    "Precomputed: No initial nucleon positions available");
		    }
		    particle.Position() =
			initial_positions[n_outgoing_hadrons_seen % initial_positions.size()];

		    event.Hadrons().push_back(particle);
		    ++n_outgoing_hadrons_seen;
		}
	    }
	}

	if(n_leptons_seen < 2) {
      	   throw std::runtime_error("Precomputed: Expected at least two leptons");
   	}

   	if(n_initial_hadrons_seen == 0) {
      	   throw std::runtime_error("Precomputed: No initial hadrons found");
   	}

	auto init_nuc = event.CurrentNucleus()->InitParticle();
	std::vector<Particle> init_parts;
	for(const auto &nucleon : event.Hadrons()) {
    	if(nucleon.Status() == ParticleStatus::initial_state) { init_parts.push_back(nucleon); }
	}

	if(init_parts.empty()) {
  	  throw std::runtime_error("Precomputed: No initial nucleons found in event");
	}

// For the event history, use one representative primary position.
// For a two-body MEC event, use the center of the sampled initial pair.
	ThreeVector primary_pos;
	for(const auto &part : init_parts) {
    	primary_pos += part.Position();
	}	
	primary_pos /= static_cast<double>(init_parts.size());

	event.History().AddVertex(primary_pos, {init_nuc}, init_parts,
                          EventHistory::StatusCode::target);

    // Setup beam in history
    auto init_lep = event.Leptons()[0];
    auto init_beam = init_lep;
    init_beam.Status() = ParticleStatus::beam;
    event.History().AddVertex({}, {init_beam}, {init_lep}, EventHistory::StatusCode::beam);

    // TODO: Figure out how to best handle tracking this with the cascade and decays
    std::vector<Particle> primary_out, propagating;
    for(const auto &part : event.Particles()) {
        if(part.IsFinal()) primary_out.push_back(part);
        if(part.IsPropagating()) {
            primary_out.push_back(part);
            propagating.push_back(part);
        }
    }
    init_parts.push_back(init_lep);
    event.History().AddVertex(init_parts[0].Position(), init_parts, primary_out,
                              EventHistory::StatusCode::primary);
}

std::unique_ptr<EventReader> NuHepMCReader::Construct(const std::string &filename,
                                                      const std::shared_ptr<Nucleus> &nuc) {
    return std::make_unique<NuHepMCReader>(filename, nuc);
}

NuHepMCReader::NuHepMCReader(const std::string &filename, std::shared_ptr<Nucleus> nuc)
    : m_nucleus{nuc} {
    spdlog::debug("[NuHepMCReader]: Initializing Reader");
    m_reader = std::make_unique<NuHepMC::Reader>(filename);
    if(!m_reader) {
        auto err_msg = fmt::format(
            "Precomputed: Failed to instantiate HepMC3::Reader from file: {}", filename);
        throw std::runtime_error(err_msg);
    }

    HepMC3::GenEvent evt;
    m_reader->read_event(evt);
    if(m_reader->failed()) {
        auto err_msg =
            fmt::format("Precomputed: Failed to read the first event from file: {}", filename);
        throw std::runtime_error(err_msg);
    }

    auto in_run_info = evt.run_info();
    if(!in_run_info) {
        // TODO: Determine to warn about missing info and not to trust header or fail
        auto err_msg = fmt::format("Precomputed: Can't read run info after first event.");
        throw std::runtime_error(err_msg);
    }

    // Read in needed header information
    m_vertex_status = NuHepMC::GR9::ReadVertexStatusIdDefinitions(in_run_info);
    m_particle_status = NuHepMC::GR10::ReadParticleStatusIdDefinitions(in_run_info);

    // TODO: Pass this to the NuHepMC Writer to ensure consistency / convert all to Achilles
    // definitions?
    // TODO: Setup Achilles components and output writer
    // TODO: Add Achilles information to the header

    // Re-open the file to ensure all events are read
    m_reader = std::make_unique<NuHepMC::Reader>(filename);
    spdlog::debug("[NuHepMCReader]: Finished initialization");
}

std::optional<achilles::Event> NuHepMCReader::Next() {
    spdlog::debug("[NuHepMCReader]: Reading next event...");
    HepMC3::GenEvent evt;
    m_reader->read_event(evt);
    if(m_reader->failed()) { return std::nullopt; }

    // TODO: Get to work for a general file (i.e. need to test on NEUT, GENIE, etc.)
    auto beampt = NuHepMC::Event::GetBeamParticle(evt);
    auto tgtpt = NuHepMC::Event::GetTargetParticle(evt);

    if(!beampt || !tgtpt) {
        spdlog::error("Precomputed: Invalid event found. Aborting");
        return std::nullopt;
    }

    // TODO: Make it work on a list of nuclei and skip instead of aborting if not available
    if(tgtpt->pdg_id() != static_cast<int>(m_nucleus->ID())) {
        spdlog::error("Precomputed: Invalid nucleus. Aborting");
        return std::nullopt;
    }

    auto primary_vtx = NuHepMC::Event::GetPrimaryVertex(evt);
    auto incoming = primary_vtx->particles_in();
    auto outgoing = primary_vtx->particles_out();

    // Convert from NuHepMC to internal event and return it
    Event event(m_nucleus, {}, 1);
    // TODO: Extend to propagate multiweights through and get the right statistics and weights
    // TODO: Handle the units to convert correctly
    event.Weight() = evt.weight();

    // Convert particles to internal ones
    // TODO: Handle multiple positions?
    ThreeVector initial_pos;
    for(auto &inpart : incoming) {
        Particle part;
        try {
            part = FromHepMC3(inpart);
        } catch(AchillesParseError &e) {
            spdlog::warn("[NuHepMCReader]: {}", e.what());
            continue;
        }
        if(part.Info().IsLepton()) event.Leptons().push_back(part);
        if(part.Info().IsHadron()) {
            event.Display();
            if(part.ID() == PID::proton()) {
                auto protons = event.Protons(ParticleStatus::background);
                auto sampled_protons = Random::Instance().Sample(1, protons);
                auto sampled_proton = sampled_protons[0];
                sampled_proton.get().Momentum() = part.Momentum();
                sampled_proton.get().Status() = ParticleStatus::initial_state;
                std::cout << sampled_proton.get() << std::endl;
                initial_pos = sampled_proton.get().Position();
            }

            if(part.ID() == PID::neutron()) {
                auto neutrons = event.Neutrons(ParticleStatus::background);
                auto sampled_neutrons = Random::Instance().Sample(1, neutrons);
                auto sampled_neutron = sampled_neutrons[0];
                sampled_neutron.get().Momentum() = part.Momentum();
                sampled_neutron.get().Status() = ParticleStatus::initial_state;
                std::cout << sampled_neutron.get() << std::endl;
                initial_pos = sampled_neutron.get().Position();
            }
        }
    }
    for(auto &outpart : outgoing) {
        Particle out;
        try {
            out = FromHepMC3(outpart);
        } catch(AchillesParseError &e) { spdlog::warn("[NuHepMCReader]: {}", e.what()); }
        if(out.Info().IsLepton()) event.Leptons().push_back(out);
        if(out.Info().IsHadron()) {
            out.Status() = ParticleStatus::propagating;
            out.Position() = initial_pos;
            event.Hadrons().push_back(out);
        }
    }
    event.Display();

    // Setup the event history
    auto init_nuc = event.CurrentNucleus()->InitParticle();
    std::vector<Particle> init_parts;
    for(const auto &nucleon : event.Hadrons()) {
        if(nucleon.Status() == ParticleStatus::initial_state) { init_parts.push_back(nucleon); }
    }
    event.History().AddVertex(init_parts[0].Position(), {init_nuc}, init_parts,
                              EventHistory::StatusCode::target);
    // Setup beam in history
    auto init_lep = event.Leptons()[0];
    auto init_beam = init_lep;
    init_beam.Status() = ParticleStatus::beam;
    event.History().AddVertex({}, {init_beam}, {init_lep}, EventHistory::StatusCode::beam);

    // TODO: Figure out how to best handle tracking this with the cascade and decays
    std::vector<Particle> primary_out, propagating;
    for(const auto &part : event.Particles()) {
        if(part.IsFinal()) primary_out.push_back(part);
        if(part.IsPropagating()) {
            primary_out.push_back(part);
            propagating.push_back(part);
        }
    }
    init_parts.push_back(init_lep);
    event.History().AddVertex(init_parts[0].Position(), init_parts, primary_out,
                              EventHistory::StatusCode::primary);

    spdlog::debug("[NuHepMCReader]: Successfully read event.");
    return event;
}

achilles::Particle NuHepMCReader::FromHepMC3(const HepMC3::ConstGenParticlePtr &part) const {
    auto momentum = part->momentum();
    FourVector mom(momentum.e(), momentum.px(), momentum.py(), momentum.pz());
    if(part->status() > 4 || part->status() < 0)
        throw AchillesParseError(
            fmt::format("Invalid particle status {}, skipping particle", part->status()));
    return Particle(part->pid(), mom, {}, static_cast<ParticleStatus>(part->status()));
}

RunCascade::RunCascade(const std::string &config_file) {
    auto config = Settings(config_file, "data/Rules_precomputed.yml");

    // Setup random number generator
    auto seed = static_cast<unsigned int>(
        std::chrono::high_resolution_clock::now().time_since_epoch().count());
    if(config.Exists("Options/Initialize/Seed"))
        if(config.GetAs<int>("Options/Initialize/Seed") > 0)
            seed = config.GetAs<unsigned int>("Options/Initialize/Seed");
    spdlog::trace("Seeding generator with: {}", seed);
    Random::Instance().Seed(seed);

    // Setup nuclei (currently only supports a single nucleus)
    auto nuclei = config.GetAs<std::vector<std::shared_ptr<Nucleus>>>("Nuclei");
    if(nuclei.size() == 0)
        throw std::runtime_error("Precomputed: Must have at least one nucleus defined");
    nucleus = nuclei[0];

    // Initialize Cascade parameters
    cascade = config.GetAs<Cascade>("Cascade");

    // Setup outputs
    bool zipped =
        config.Exists("Main/Output/Zipped") ? config.GetAs<bool>("Main/Output/Zipped") : true;
    auto format = config.GetAs<std::string>("Main/Output/Format");
    auto name = config.GetAs<std::string>("Main/Output/Name");
    spdlog::trace("Outputing as {} format", format);
    if(format == "Achilles") {
        writer = std::make_unique<AchillesWriter>(name, zipped);
#ifdef ACHILLES_ENABLE_HEPMC3
    } else if(format == "HepMC3") {
        writer = std::make_unique<HepMC3Writer>(name, zipped);
    } else if(format == "NuHepMC") {
        writer = std::make_unique<NuHepMCWriter>(name, zipped);
#endif
    } else {
        std::string msg = fmt::format("Achilles: Invalid output format requested {}", format);
        throw std::runtime_error(msg);
    }
    std::vector<ProcessGroup> dummy;
    writer->WriteHeader(config_file, dummy);

    // Setup input
    format = config.GetAs<std::string>("Main/Input/Format");
    auto filename = config.GetAs<std::string>("Main/Input/Name");
    reader = EventReaderFactory::Initialize(format, filename, nucleus);
}

void RunCascade::RunAll() {
    auto run = [&](Event &event) { Run(event); };
    EventPipeline pipeline(run);
    RunPipeline(reader, pipeline);
}

void RunCascade::Run(Event &event) {
    cascade.Evolve(event, nucleus.get());
    std::vector<Particle> init_part, final_part;
    for(const auto &part : event.Particles()) {
        if(part.IsFinal()) final_part.push_back(part);
        if(part.IsInitial()) init_part.push_back(part);
    }
    event.History().UpdateStatuses(event.Hadrons());
    writer->Write(event);
}
