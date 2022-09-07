#include "Achilles/ProcessInfo.hh"

extern "C" {
    void ProcessIDs(achilles::Process_Info *info, long *ids, size_t &len) {
        len = info -> ids.size();
        ids = new long[info -> ids.size()];
        for(size_t i = 0; i < info -> ids.size(); ++i) {
            ids[i] = info -> ids[i].AsInt();
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
        info -> state = {part_in, part_out};
    }
}
