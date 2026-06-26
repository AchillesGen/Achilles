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
	size_t size() const override { return lbeam->size() + hbeam->size() + finalstate->size(); }
    size_t NDims() const override { return lbeam->NDims() + hbeam->NDims() + finalstate->NDims(); }
    void SetLeptonBeam(Mapper_sptr<FourVector> _lbeam) { lbeam = _lbeam; }
    void SetHadronBeam(Mapper_sptr<FourVector> _hbeam) { hbeam = _hbeam; }
    void SetFinalState(Mapper_ptr<FourVector> _finalstate) { finalstate = std::move(_finalstate); }
	/// Start of Lepton list
	size_t LeptonIdx() const { return 0; }
	/// Start of Hadron list
	size_t HadronIdx() const {
		try {
			return lbeam->size()+LeptonIdx();
		} catch(...) {
			// Assume lbeam has a size of 1 if not initialized
			return 1+LeptonIdx();
		}
	}
	/// Start of Final State
	size_t FinalStateIdx() const {
		try {
			return hbeam->size()+HadronIdx();
		} catch(...) {
			return 1+HadronIdx();
		}
	}

    PSMapper(size_t _nleptons, size_t _nhadrons, size_t _nspectators)
        : nleptons{std::move(_nleptons)}, nhadrons{std::move(_nhadrons)},
          nspectators{_nspectators} {}

  private:
    size_t nleptons, nhadrons, nspectators;
    Mapper_sptr<FourVector> lbeam{}, hbeam{};
    Mapper_ptr<FourVector> finalstate{};
};

} // namespace achilles

#endif
