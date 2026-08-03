#ifndef BEAM_MAPPER_HH
#define BEAM_MAPPER_HH

#include "Achilles/Mapper.hh"

namespace achilles {

class Beam;
class Event;

class BeamMapper : public Mapper<Event> {
  public:
    BeamMapper(std::shared_ptr<Beam> beam)
        : m_beam{std::move(beam)} {}

    void GeneratePoint(Event&, const std::vector<double> &) override;
    double GenerateWeight(const Event&, std::vector<double> &) override;
    size_t NDims() const override;

  private:
    std::shared_ptr<Beam> m_beam;
};

} // namespace achilles

#endif
