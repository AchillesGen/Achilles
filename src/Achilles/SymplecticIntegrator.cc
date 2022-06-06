#include "Achilles/SymplecticIntegrator.hh"

using SI = achilles::SymplecticIntegrator;

void SI::Initialize(const ThreeVector &q, const ThreeVector &p) {
    m_state = PSState(q, p);
}

void SI::HamiltonianA(double time_step) {
    m_state.p -= time_step*m_dHdr(m_state.q, m_state.y, m_pot);
    m_state.x += time_step*m_dHdp(m_state.q, m_state.y, m_pot);
}

void SI::HamiltonianB(double time_step) {
    m_state.q += time_step*m_dHdp(m_state.x, m_state.p, m_pot);
    m_state.y -= time_step*m_dHdr(m_state.x, m_state.p, m_pot);
}

void SI::Coupling(double time_step) {
    const double comega = cos(2*m_omega*time_step);
    const double somega = sin(2*m_omega*time_step);
    const auto qsum = m_state.q + m_state.x;
    const auto psum = m_state.p + m_state.y;
    const auto qdiff = m_state.q - m_state.x;
    const auto pdiff = m_state.p - m_state.y;

    m_state.q = (qsum + comega*qdiff + somega*pdiff)/2;
    m_state.p = (psum - somega*qdiff + comega*pdiff)/2;
    m_state.x = (qsum - comega*qdiff - somega*pdiff)/2;
    m_state.y = (psum + somega*qdiff - comega*pdiff)/2;
}
