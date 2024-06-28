#include "Achilles/ProcessInfo.hh"

extern "C" {
void ProcessIDs(achilles::ProcessInfo *info, long *ids, size_t &len) {
    auto pids = info->Ids();
    len = pids.size();
    ids = pids.data();
}
size_t ProcessMultiplicity(achilles::ProcessInfo *info) {
    return info->Multiplicity();
}
void ProcessMasses(achilles::ProcessInfo *info, double *masses, size_t &len) {
    masses = info->Masses().data();
    len = info->Masses().size();
}
void ProcessAddState(achilles::ProcessInfo *info, long *initial, long *final, size_t in,
                     size_t out) {
    std::vector<achilles::PID> part_in, part_out;
    for(size_t i = 0; i < in; ++i) part_in.push_back(achilles::PID(initial[i]));
    for(size_t i = 0; i < out; ++i) part_out.push_back(achilles::PID(final[i]));
    info->m_hadronic = {part_in, part_out};
}
}
