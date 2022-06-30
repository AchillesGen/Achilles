#include "nuchic/ProcessInfo.hh"

extern "C" {
    const char* ProcessModel(nuchic::Process_Info *info) { return info -> m_model.c_str(); }
    void ProcessIDs(nuchic::Process_Info *info, long *ids, size_t &len) {
        len = info -> m_ids.size();
        ids = new long[info -> m_ids.size()];
        for(size_t i = 0; i < info -> m_ids.size(); ++i) {
            ids[i] = info -> m_ids[i].AsInt();
        }
    }
    size_t ProcessMultiplicity(nuchic::Process_Info *info) { return info -> Multiplicity(); }
    void ProcessMasses(nuchic::Process_Info *info, double *masses, size_t &len) {
        masses = info -> Masses().data();
        len = info -> Masses().size();
    }
    void ProcessAddState(nuchic::Process_Info *info, long *initial, long *final, size_t in, size_t out) {
        std::vector<nuchic::PID> part_in, part_out;
        for(size_t i = 0; i < in; ++i) part_in.push_back(nuchic::PID(initial[i]));
        for(size_t i = 0; i < out; ++i) part_out.push_back(nuchic::PID(final[i]));
        info -> m_states[part_in] = part_out;
    }
}
