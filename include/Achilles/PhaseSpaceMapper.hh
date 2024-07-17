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
    void SetLeptonBeam(Mapper_sptr<FourVector> _lbeam) { lbeam = _lbeam; }
    void SetHadronBeam(Mapper_sptr<FourVector> _hbeam) { hbeam = _hbeam; }
    void SetFinalState(Mapper_ptr<FourVector> final) { main = std::move(final); }

    YAML::Node ToYAML() const override {
        YAML::Node node;
        node["Name"] = "PSMapper";
        node["nlep"] = nleptons;
        node["nhad"] = nhadrons;
        node["BeamMapper"] = lbeam->ToYAML();
        node["HadronMapper"] = hbeam->ToYAML();
        node["FSMapper"] = main->ToYAML();
        return node;
    }
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
