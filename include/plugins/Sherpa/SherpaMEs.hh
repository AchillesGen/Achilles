#ifndef achilles__plugins__SherpaMEs_hh
#define achilles__plugins__SherpaMEs_hh

#include "Achilles/Achilles.hh"
#include "Achilles/ParticleInfo.hh"

#include <complex>
#include <array>
#include <vector>
#include <string>

namespace ATOOLS { class Cluster_Amplitude; }
namespace PHASIC { 
    class Process_Base;
    class Channels;
}
namespace SHERPA { class Sherpa; }

namespace achilles {

struct Process_Info;
struct FormFactorInfo;

class SherpaMEs {
private:

  SHERPA::Sherpa *p_sherpa{};

  void addParameter(std::vector<char*>& argv,const std::string& val) const;
  int SherpaVerbosity(int loglevel) const;

  PHASIC::Process_Base *getProcess(ATOOLS::Cluster_Amplitude* const ampl);

public:

  using LeptonCurrents = std::map<int, std::vector<std::vector<std::complex<double>>>>;

  SherpaMEs() = default;

  MOCK ~SherpaMEs();

  bool Initialize(const std::vector<std::string> &args);
  bool InitializeProcess(const Process_Info &info);

  std::vector<std::unique_ptr<PHASIC::Channels>> GenerateChannels(const std::vector<long> &fl) const;
  std::map<size_t, long> MomentumMap(const std::vector<long> &fl) const;

  MOCK LeptonCurrents Calc
  (const std::vector<int> &fl,
   const std::vector<std::array<double, 4> > &p,
   const double &mu2) const;

  MOCK std::vector<FormFactorInfo> FormFactors(int, int) const;
  double Coupling(const std::string&) const;
  void RegisterParticles() const;

};// end of class SherpaMEs

}// end namespace achilles

#endif
