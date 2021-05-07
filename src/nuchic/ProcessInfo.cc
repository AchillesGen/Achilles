#include "nuchic/ProcessInfo.hh"
#include "nuchic/Beams.hh"

void nuchic::Process_Info::AddBeam(const nuchic::Beam &beam) {
    m_ids.insert(m_ids.begin(), nuchic::PID::dummyHadron());
    m_ids.insert(m_ids.begin(), *beam.BeamIDs().begin());
}
