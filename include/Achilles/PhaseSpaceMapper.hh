#ifndef PHASESPACE_MAPPER_HH
#define PHASESPACE_MAPPER_HH

#include "Achilles/Mapper.hh"
#include "Achilles/Event.hh"

namespace achilles {

class FourVector;

class PSBuilder;

class PSMapper : public Mapper<Event> {
  public:
    friend class PSBuilder;

    void GeneratePoint(Event&, const std::vector<double> &) override;
    double GenerateWeight(const Event&, std::vector<double> &) override;
    size_t NDims() const override { return lbeam->NDims() + hbeam->NDims() + finalstate->NDims(); }
    void SetLeptonBeam(Mapper_sptr<Event> _lbeam) { lbeam = _lbeam; }
    void SetHadronBeam(Mapper_sptr<Event> _hbeam) { hbeam = _hbeam; }
    void SetFinalState(Mapper_ptr<Event> _finalstate) { finalstate = std::move(_finalstate); }

    PSMapper(size_t _nleptons, size_t _nhadrons, size_t _nspectators)
        : nleptons{std::move(_nleptons)}, nhadrons{std::move(_nhadrons)},
          nspectators{_nspectators} {}

  private:
    size_t nleptons, nhadrons, nspectators;
    Mapper_sptr<Event> lbeam{}, hbeam{};
    Mapper_ptr<Event> finalstate{};
};

} // namespace achilles

#endif