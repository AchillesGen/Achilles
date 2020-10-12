#include "nuchic/Beams.hh"


nuchic::Beam::Beam(BeamMap beams) : m_beams{std::move(beams)} {
    n_vars = 0;
    for(const auto& beam : m_beams) {
        if(beam.second -> NVariables() > n_vars)
            n_vars = beam.second -> NVariables();
        if(m_pids.find(beam.first) != m_pids.end())
            throw std::logic_error(fmt::format("Multiple beams exist for PID: {}", int(beam.first)));
        m_pids.insert(beam.first);
    }
}
