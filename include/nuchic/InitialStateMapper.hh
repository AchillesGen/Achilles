#ifndef INITIALSTATE_MAPPER_HH
#define INITIALSTATE_MAPPER_HH

#include "nuchic/Mapper.hh"

namespace nuchic {

class FourVector;
class Beam;

class InitialMapper : public Mapper<FourVector> {
    public:
        InitialMapper(std::shared_ptr<Beam> beam) : m_beam{std::move(beam)} {}

        void GeneratePoint(std::vector<FourVector>&, const std::vector<double>&) const override;
        double GenerateWeight(const std::vector<FourVector>&, std::vector<double>&) const override;
        size_t NDims() const override;

    private:
        std::shared_ptr<Beam> m_beam;
};

}

#endif
