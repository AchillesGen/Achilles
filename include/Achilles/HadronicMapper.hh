#ifndef HADRONIC_MAPPER_HH
#define HADRONIC_MAPPER_HH

#include <cmath>

#include "Achilles/Mapper.hh"
#include "Achilles/PhaseSpaceFactory.hh"

namespace achilles {

class FourVector;

class HadronicBeamMapper : public Mapper<FourVector> {
  public:
    HadronicBeamMapper(size_t idx, std::string name)
        : m_idx{std::move(idx)}, m_name{std::move(name)} {}

    void GeneratePoint(std::vector<FourVector> &, const std::vector<double> &) override = 0;
    double GenerateWeight(const std::vector<FourVector> &, std::vector<double> &) override = 0;
    size_t NDims() const override = 0;
    YAML::Node ToYAML() const override {
        YAML::Node result;
        result["Name"] = m_name;
        result["Idx"] = m_idx;
        return result;
    }
    static std::string Name() { return "Hadronic Initial State"; }

  protected:
    size_t HadronIdx() const { return m_idx; }

  private:
    size_t m_idx;
    std::string m_name;
};

class QESpectralMapper : public HadronicBeamMapper,
                         RegistrablePS<HadronicBeamMapper, QESpectralMapper, size_t> {
  public:
    QESpectralMapper(size_t idx) : HadronicBeamMapper(idx, Name()) {}
    static std::string Name() { return "QESpectral"; }
    static std::unique_ptr<HadronicBeamMapper> Construct(const size_t &idx) {
        return std::make_unique<QESpectralMapper>(idx);
    }

    void GeneratePoint(std::vector<FourVector> &, const std::vector<double> &) override;
    double GenerateWeight(const std::vector<FourVector> &, std::vector<double> &) override;
    size_t NDims() const override { return 4; }

  private:
    // static constexpr double dCos = 2;
    static constexpr double dPhi = 2 * M_PI;
    // static constexpr double dp = 800;
    // static constexpr double dE = 400;
};

class CoherentMapper : public HadronicBeamMapper,
                       RegistrablePS<HadronicBeamMapper, CoherentMapper, size_t> {
  public:
    CoherentMapper(size_t idx) : HadronicBeamMapper(idx, Name()) {}
    static std::string Name() { return "Coherent"; }
    static std::unique_ptr<HadronicBeamMapper> Construct(const size_t &idx) {
        return std::make_unique<CoherentMapper>(idx);
    }

    void GeneratePoint(std::vector<FourVector> &, const std::vector<double> &) override;
    double GenerateWeight(const std::vector<FourVector> &, std::vector<double> &) override;
    size_t NDims() const override { return 0; }
};


class IntfSpectralMapper : public HadronicBeamMapper,
                         RegistrablePS<HadronicBeamMapper, IntfSpectralMapper, size_t> {
  public:
    IntfSpectralMapper(size_t idx) : HadronicBeamMapper(idx, Name()) {}
    static std::string Name() { return "IntfSpectral"; }
    static std::unique_ptr<HadronicBeamMapper> Construct(const size_t &idx) {
        return std::make_unique<IntfSpectralMapper>(idx);
    }

    void GeneratePoint(std::vector<FourVector> &, const std::vector<double> &) override;
    double GenerateWeight(const std::vector<FourVector> &, std::vector<double> &) override;
    size_t NDims() const override { return 7; }

  private:
    static constexpr double dCos2 = 2;
    static constexpr double dPhi = 2 * M_PI;
    static constexpr double dp2 = 400; // reduced b/c MF SF has no strength at high P
    // static constexpr double dE = 400;
};

class MecSpectralMapper : public HadronicBeamMapper,
                         RegistrablePS<HadronicBeamMapper, MecSpectralMapper, size_t> {
  public:
    MecSpectralMapper(size_t idx) : HadronicBeamMapper(idx, Name()) {}
    static std::string Name() { return "MecSpectral"; }
    static std::unique_ptr<HadronicBeamMapper> Construct(const size_t &idx) {
        return std::make_unique<MecSpectralMapper>(idx);
    }

    void GeneratePoint(std::vector<FourVector> &, const std::vector<double> &) override;
    double GenerateWeight(const std::vector<FourVector> &, std::vector<double> &) override;
    size_t NDims() const override { return 8; }

  private:
    static constexpr double dCos2 = 2;
    static constexpr double dPhi = 2 * M_PI;
    static constexpr double dp2 = 800; 
    static constexpr double dE2 = 80;
};



} // namespace achilles

#endif
