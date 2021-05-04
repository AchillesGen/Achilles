#ifndef nuchic__plugins__SherpaMEs_hh
#define nuchic__plugins__SherpaMEs_hh

#include "nuchic/ParticleInfo.hh"

#include <vector>
#include <string>

namespace ATOOLS { class Cluster_Amplitude; }
namespace PHASIC { class Process_Base; }
namespace SHERPA { class Sherpa; }

namespace nuchic {

struct Process_Info {
  std::vector<nuchic::PID> m_ids;
  inline Process_Info(const std::vector<nuchic::PID> &ids=
		      std::vector<nuchic::PID>()):
    m_ids(ids) {}
};// end of struct Process_Info 

class SherpaMEs {
private:

  SHERPA::Sherpa *p_sherpa;

  void addParameter(std::vector<char*>& argv,const std::string& val) const;
  int SherpaVerbosity(int loglevel) const;

  PHASIC::Process_Base *getProcess(ATOOLS::Cluster_Amplitude* const ampl);

public:

  SherpaMEs(): p_sherpa(NULL) {}

  ~SherpaMEs();

  bool Initialize(const std::vector<std::string> &args);

  bool InitializeProcess(const Process_Info &info);

};// end of class SherpaMEs

}// end namespace nuchic

#endif
