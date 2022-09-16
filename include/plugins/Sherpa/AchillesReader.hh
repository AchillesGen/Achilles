#pragma once
#include "PHASIC++/Main/Event_Reader.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/Message.H"

using namespace PHASIC;
using namespace ATOOLS;

namespace achilles {

  class Achilles_Reader: public Event_Reader {
  private:

    Cluster_Amplitude *p_ampl;
    
  public:

    Achilles_Reader(const Event_Reader_Key &key):
      Event_Reader(key), p_ampl(NULL) {}

    ~Achilles_Reader() {}

    Cluster_Amplitude *ReadEvent() {
        return p_ampl;
    }

    void SetAmpl(Cluster_Amplitude *ampl) {
        p_ampl = ampl;
    }
    
  };// end of class Achilles_Reader
  
}// end of namespace achilles

using namespace achilles;

DECLARE_GETTER(Achilles_Reader,"Achilles",Event_Reader,Event_Reader_Key);

Event_Reader *ATOOLS::Getter<Event_Reader,Event_Reader_Key,Achilles_Reader>::
operator()(const Event_Reader_Key &args) const
{
  return new Achilles_Reader(args);
}

void ATOOLS::Getter<Event_Reader,Event_Reader_Key,Achilles_Reader>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"Achilles reader";
}
