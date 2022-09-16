#ifndef achilles__plugins__SherpaMEs_hh
#define achilles__plugins__SherpaMEs_hh

#include "Achilles/Achilles.hh"
#include "Achilles/ParticleInfo.hh"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Woverloaded-virtual"
#pragma GCC diagnostic ignored "-Wignored-qualifiers"
#pragma GCC diagnostic ignored "-Wuseless-cast"
#include "COMIX/Main/Single_Process.H"
#include "plugins/Sherpa/AchillesReader.hh"
#pragma GCC diagnostic pop

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
  COMIX::Single_Process *singleProcess;

  achilles::Achilles_Reader *reader;

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
   const double &mu2);
  void FillAmplitudes(std::vector<Spin_Amplitudes> &amps);

  MOCK std::vector<FormFactorInfo> FormFactors(int, int) const;
  double Coupling(const std::string&) const;
  void RegisterParticles() const;

};// end of class SherpaMEs

}// end namespace achilles

#endif
