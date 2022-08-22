#include "Achilles/ProcessInfo.hh"

extern "C" {
    const char* ProcessModel(achilles::Process_Info *info) { return info -> m_model.c_str(); }
    void ProcessIDs(achilles::Process_Info *info, long *ids, size_t &len) {
        len = info -> m_ids.size();
        ids = new long[info -> m_ids.size()];
        for(size_t i = 0; i < info -> m_ids.size(); ++i) {
            ids[i] = info -> m_ids[i].AsInt();
        }
    }
    size_t ProcessMultiplicity(achilles::Process_Info *info) { return info -> Multiplicity(); }
    void ProcessMasses(achilles::Process_Info *info, double *masses, size_t &len) {
        masses = info -> Masses().data();
        len = info -> Masses().size();
    }
    void ProcessAddState(achilles::Process_Info *info, long *initial, long *final, size_t in, size_t out) {
        std::vector<achilles::PID> part_in, part_out;
        for(size_t i = 0; i < in; ++i) part_in.push_back(achilles::PID(initial[i]));
        for(size_t i = 0; i < out; ++i) part_out.push_back(achilles::PID(final[i]));
        info -> m_states[part_in] = part_out;
    }
}
