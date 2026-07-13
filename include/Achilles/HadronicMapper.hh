#ifndef HADRONIC_MAPPER_HH
#define HADRONIC_MAPPER_HH

#include <cmath>

#include "Achilles/Factory.hh"
#include "Achilles/Mapper.hh"
#include "Achilles/ProcessInfo.hh"

namespace achilles {

class FourVector;
class Event;

class HadronicBeamMapper : public Mapper<Event> {
  public:
    HadronicBeamMapper(const ProcessInfo &info, size_t idx) : m_info{info}, m_idx{idx} {}

    void GeneratePoint(Event&, const std::vector<double> &) override = 0;
    double GenerateWeight(const Event&, std::vector<double> &) override = 0;
    size_t NDims() const override = 0;
    static std::string Name() { return "Hadronic Initial State"; }

  protected:
    size_t HadronIdx() const { return m_idx; }
    ProcessInfo m_info;

  private:
    size_t m_idx;
};

class QESpectralMapper
    : public HadronicBeamMapper,
      Registrable<HadronicBeamMapper, QESpectralMapper, const ProcessInfo &, size_t> {
  public:
    QESpectralMapper(const ProcessInfo &info, size_t idx);
    static std::string Name() { return "OneBodySpectral"; }
    static std::unique_ptr<HadronicBeamMapper> Construct(const ProcessInfo &info, size_t idx) {
        return std::make_unique<QESpectralMapper>(info, idx);
    }

    void GeneratePoint(Event&, const std::vector<double> &) override;
    double GenerateWeight(const Event&, std::vector<double> &) override;
    size_t NDims() const override { return 4; }

  private:
    // static constexpr double dCos = 2;
    static constexpr double dPhi = 2 * M_PI;
    // static constexpr double dp = 800;
    // static constexpr double dE = 400;
};

class CoherentMapper
    : public HadronicBeamMapper,
      Registrable<HadronicBeamMapper, CoherentMapper, const ProcessInfo &, size_t> {
  public:
    CoherentMapper(const ProcessInfo &info, size_t idx);
    static std::string Name() { return "Coherent"; }
    static std::unique_ptr<HadronicBeamMapper> Construct(const ProcessInfo &info, size_t idx) {
        return std::make_unique<CoherentMapper>(info, idx);
    }

    void GeneratePoint(Event&, const std::vector<double> &) override;
    double GenerateWeight(const Event&, std::vector<double> &) override;
    size_t NDims() const override { return 0; }

  private:
    double m_mass;
};

class IntfSpectralMapper
    : public HadronicBeamMapper,
      Registrable<HadronicBeamMapper, IntfSpectralMapper, const ProcessInfo &, size_t> {
  public:
    IntfSpectralMapper(const ProcessInfo &info, size_t idx);
    static std::string Name() { return "IntfSpectral"; }
    static std::unique_ptr<HadronicBeamMapper> Construct(const ProcessInfo &info, size_t idx) {
        return std::make_unique<IntfSpectralMapper>(info, idx);
    }

    void GeneratePoint(Event&, const std::vector<double> &) override;
    double GenerateWeight(const Event&, std::vector<double> &) override;
    size_t NDims() const override { return 7; }

  private:
    static constexpr double dCos2 = 2;
    static constexpr double dPhi = 2 * M_PI;
    static constexpr double dp2 = 400; // reduced b/c MF SF has no strength at high P
    // static constexpr double dE = 400;
};

class DISSingleNucleonMapper: public HadronicBeamMapper,
		Registrable<HadronicBeamMapper,DISSingleNucleonMapper,const ProcessInfo&,size_t> {
  public:
	DISSingleNucleonMapper(const ProcessInfo&,size_t);
	static std::string Name() { return "DISSingleNucleon"; }
	static std::unique_ptr<HadronicBeamMapper> Construct(const ProcessInfo& info,size_t idx) {
		return std::make_unique<DISSingleNucleonMapper>(info,idx);
	}
	void GeneratePoint(Event&,const std::vector<double>&) override;
	double GenerateWeight(const Event&,std::vector<double>&) override;

	/// 1 for [x]
	size_t NDims() const override { return 1; }
  private:
	FourVector pHadron;
};

class DISNucleusMapper: public HadronicBeamMapper,
		Registrable<HadronicBeamMapper,DISNucleusMapper,const ProcessInfo&,size_t> {
  public:
	DISNucleusMapper(const ProcessInfo&,size_t);
	static std::string Name() { return "DISNucleus"; }
	static std::unique_ptr<HadronicBeamMapper> Construct(const ProcessInfo& info,size_t idx) {
		return std::make_unique<DISSingleNucleonMapper>(info,idx);
	}
	void GeneratePoint(Event&,const std::vector<double>&) override;
	double GenerateWeight(const Event&,std::vector<double>&) override;

	/// 2 for [pHadron,pQuark]
	size_t size() const override { return 2; }
	/// 5 for [|pHadron|,cosTheta,phi,energy,x]
	size_t NDims() const override { return 5; }
  private:
	static constexpr double dPhi=2*M_PI;
};

} // namespace achilles

#endif