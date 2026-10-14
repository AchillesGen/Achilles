#ifndef PHASESPACE_MAPPER_HH
#define PHASESPACE_MAPPER_HH

#include "Achilles/Mapper.hh"

namespace achilles {

class FourVector;

class PSBuilder;

class PSMapper : public Mapper<FourVector> {
  public:
    friend class PSBuilder;
    static PSBuilder build(size_t, size_t, size_t);

    void GeneratePoint(std::vector<FourVector> &, const std::vector<double> &) override;
    double GenerateWeight(const std::vector<FourVector> &, std::vector<double> &) override;
    size_t NDims() const override { return lbeam->NDims() + hbeam->NDims() + main->NDims(); }
    // The random numbers are ordered [hbeam][lbeam][main], so a dimension the
    // lepton beam fixes sits after the hadron beam's.
    std::optional<size_t> FixedDim() const override {
        const auto dim = lbeam->FixedDim();
        if(!dim) return std::nullopt;
        return hbeam->NDims() + *dim;
    }
    double FixedDensity() const override { return lbeam->FixedDensity(); }
    void SetLeptonBeam(Mapper_sptr<FourVector> _lbeam) { lbeam = _lbeam; }
    void SetHadronBeam(Mapper_sptr<FourVector> _hbeam) { hbeam = _hbeam; }
    void SetFinalState(Mapper_ptr<FourVector> final) { main = std::move(final); }

    PSMapper(size_t _nleptons, size_t _nhadrons, size_t _nspectators)
        : nleptons{std::move(_nleptons)}, nhadrons{std::move(_nhadrons)},
          nspectators{_nspectators} {}

  private:
    size_t nleptons, nhadrons, nspectators;
    Mapper_sptr<FourVector> lbeam{}, hbeam{};
    Mapper_ptr<FourVector> main{};
};

} // namespace achilles

#endif
