#ifndef PLUGINS_CHANNELS_HH
#define PLUGINS_CHANNELS_HH

#include "ATOOLS/Math/Vector.H"
#include "nuchic/Mapper.hh"

namespace PHASIC {

class Channels : public nuchic::Mapper<ATOOLS::Vec4D> {
    public:
        double GenerateWeight(const std::vector<ATOOLS::Vec4D>&, std::vector<double>&) const override = 0;
        void GeneratePoint(std::vector<ATOOLS::Vec4D>&, const std::vector<double>&) const override = 0;
        size_t NDims() const override = 0;
        static std::string Name() { return "Sherpa Final State"; }
};

}

#endif
