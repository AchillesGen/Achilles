#ifndef BEAM_MAPPER_HH
#define BEAM_MAPPER_HH

#include "nuchic/Mapper.hh"

namespace nuchic {

class FourVector;
class Beam;

class BeamMapper : public Mapper<FourVector> {
    public:
        BeamMapper(size_t idx, std::shared_ptr<Beam> beam) : m_idx{std::move(idx)}, m_beam{std::move(beam)} {}

        void GeneratePoint(std::vector<FourVector>&, const std::vector<double>&) const override;
        double GenerateWeight(const std::vector<FourVector>&, std::vector<double>&) const override;
        size_t NDims() const override;
        YAML::Node ToYAML() const override {
            YAML::Node result;
            result["Name"] = "Beam";
            result["Beam"] = *m_beam;
            return result;
        }

    private:
        size_t m_idx;
        std::shared_ptr<Beam> m_beam;
};

}

#endif
