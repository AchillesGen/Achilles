#ifndef achilles__plugins__SherpaInterface_hh
#define achilles__plugins__SherpaInterface_hh

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wuseless-cast"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Math/Vector.H"
#pragma GCC diagnostic pop
#include "Achilles/Achilles.hh"
#include "Achilles/Particle.hh"
#include "Achilles/ParticleInfo.hh"
#include "Achilles/EventHistory.hh"

#include <complex>
#include <array>
#include <vector>
#include <string>

namespace ATOOLS { 
    class Cluster_Amplitude;
    class Particle;
}
namespace PHASIC { 
    class Process_Base;
    class Channels;
}
namespace SHERPA { class Sherpa; }
namespace COMIX { class Single_Process; }
namespace METOOLS { class Spin_Amplitudes; }

namespace achilles {

struct Process_Info;
struct FormFactorInfo;
class Achilles_Reader;
class Event;

class SherpaInterface {
private:

  SHERPA::Sherpa *p_sherpa{};

  void addParameter(std::vector<char*>& argv,const std::string& val) const;
  int SherpaVerbosity(int loglevel) const;
  static FourVector ToAchilles(const ATOOLS::Vec4D&);
  static Particle ToAchilles(ATOOLS::Particle*, bool);
  static void AddHistoryNode(ATOOLS::Blob* blob, EventHistory &history, EventHistory::StatusCode);
  static void ToAchilles(ATOOLS::Blob_List*, EventHistory&);
  static ATOOLS::Blob_List FromAchilles(EventHistory);

  PHASIC::Process_Base *getProcess(ATOOLS::Cluster_Amplitude* const ampl);
  COMIX::Single_Process *singleProcess;
  achilles::Achilles_Reader *reader;


public:

  using LeptonCurrents = std::map<int, std::vector<std::vector<std::complex<double>>>>;

  SherpaInterface() = default;

  MOCK ~SherpaInterface();

  bool Initialize(const std::vector<std::string> &args);
  bool InitializeProcess(const Process_Info &info);

  std::vector<std::unique_ptr<PHASIC::Channels>> GenerateChannels(const std::vector<long> &fl) const;
  std::map<size_t, long> MomentumMap(const std::vector<long> &fl) const;

  MOCK LeptonCurrents Calc
  (const std::vector<int> &fl,
   const std::vector<std::array<double, 4> > &p,
   const double &mu2);
  void FillAmplitudes(std::vector<METOOLS::Spin_Amplitudes> &amps);

  MOCK std::vector<FormFactorInfo> FormFactors(int, int) const;
  double Coupling(const std::string&) const;
  void RegisterParticles() const;

  void GenerateEvent(Event&);
};// end of class SherpaInterface

}// end namespace achilles

#endif
