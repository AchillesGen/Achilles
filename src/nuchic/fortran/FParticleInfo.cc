#include "nuchic/ParticleInfo.hh"

extern "C" {

    using PInfo = nuchic::ParticleInfo;

    PInfo* CreateParticleInfo(const long int id) {
        return new PInfo(id);
    }

    void DeleteParticleInfo(PInfo *self) {
        delete self;
    }

    // Functions for ParticleInfo
    char* Name(const PInfo *self) {
        auto name = self -> Name();
        size_t len = name.size();
        char* cname = new char[len];
        strcpy(cname, self -> Name().c_str());
        return cname;
    }

    int PID(const PInfo *self) {
        return self -> IntID();
    }

    double Charge(const PInfo *self) {
        return self -> Charge();
    }

    double Spin(const PInfo *self) {
        return self -> Spin();
    }

    double Mass(const PInfo *self) {
        return self -> Mass();
    }

    double Width(const PInfo *self) {
        return self -> Width();
    }
}
