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
      Event_Reader(key), p_ampl(NULL) 
    {
      m_compute = false;
    }

    ~Achilles_Reader() {}

    Cluster_Amplitude *ReadEvent() {
        return p_ampl;
    }

    void SetAmpl(Cluster_Amplitude *ampl) {
        p_ampl = ampl;
        // std::swap<Cluster_Leg*>(p_ampl->Legs()[0], p_ampl->Legs()[1]);
    }

    Cluster_Amplitude* GetAmpl() { return p_ampl; }
    
  };// end of class Achilles_Reader
  
}// end of namespace achilles
