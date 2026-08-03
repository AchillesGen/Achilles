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
    HadronicBeamMapper(const ProcessInfo &info) : m_info{info} {}

    void GeneratePoint(Event&, const std::vector<double> &) override = 0;
    double GenerateWeight(const Event&, std::vector<double> &) override = 0;
    size_t NDims() const override = 0;
    static std::string Name() { return "Hadronic Initial State"; }

  protected:
    ProcessInfo m_info;
};

class QESpectralMapper
    : public HadronicBeamMapper,
      Registrable<HadronicBeamMapper, QESpectralMapper, const ProcessInfo &> {
  public:
    QESpectralMapper(const ProcessInfo &info);
    static std::string Name() { return "OneBodySpectral"; }
    static std::unique_ptr<HadronicBeamMapper> Construct(const ProcessInfo &info) {
        return std::make_unique<QESpectralMapper>(info);
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
      Registrable<HadronicBeamMapper, CoherentMapper, const ProcessInfo &> {
  public:
    CoherentMapper(const ProcessInfo &info);
    static std::string Name() { return "Coherent"; }
    static std::unique_ptr<HadronicBeamMapper> Construct(const ProcessInfo &info) {
        return std::make_unique<CoherentMapper>(info);
    }

    void GeneratePoint(Event&, const std::vector<double> &) override;
    double GenerateWeight(const Event&, std::vector<double> &) override;
    size_t NDims() const override { return 0; }

  private:
    double m_mass;
};

class IntfSpectralMapper
    : public HadronicBeamMapper,
      Registrable<HadronicBeamMapper, IntfSpectralMapper, const ProcessInfo &> {
  public:
    IntfSpectralMapper(const ProcessInfo &info);
    static std::string Name() { return "IntfSpectral"; }
    static std::unique_ptr<HadronicBeamMapper> Construct(const ProcessInfo &info) {
        return std::make_unique<IntfSpectralMapper>(info);
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
		Registrable<HadronicBeamMapper,DISSingleNucleonMapper,const ProcessInfo&> {
  public:
	DISSingleNucleonMapper(const ProcessInfo&);
	static std::string Name() { return "DISSingleNucleon"; }
	static std::unique_ptr<HadronicBeamMapper> Construct(const ProcessInfo& info) {
		return std::make_unique<DISSingleNucleonMapper>(info);
	}
	void GeneratePoint(Event&,const std::vector<double>&) override;
	double GenerateWeight(const Event&,std::vector<double>&) override;

	/// 1 for [x]
	size_t NDims() const override { return 1; }
  private:
	FourVector pHadron;
};

class DISNucleusMapper: public HadronicBeamMapper,
		Registrable<HadronicBeamMapper,DISNucleusMapper,const ProcessInfo&> {
  public:
	DISNucleusMapper(const ProcessInfo&);
	static std::string Name() { return "DISNucleus"; }
	static std::unique_ptr<HadronicBeamMapper> Construct(const ProcessInfo& info) {
		return std::make_unique<DISNucleusMapper>(info);
	}
	void GeneratePoint(Event&,const std::vector<double>&) override;
	double GenerateWeight(const Event&,std::vector<double>&) override;

	/// 5 for [|pHadron|,cosTheta,phi,energy,x]
	size_t NDims() const override { return 5; }
  private:
	static constexpr double dPhi=2*M_PI;
};

} // namespace achilles

#endif