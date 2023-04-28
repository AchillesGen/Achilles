#pragma once
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "PHASIC++/Main/Event_Reader.H"

using namespace PHASIC;
using namespace ATOOLS;

namespace achilles {

class Achilles_Reader : public Event_Reader {
  public:
    Achilles_Reader(const Event_Reader_Key &key) : Event_Reader(key) { m_compute = false; }

    ~Achilles_Reader() {}

    Cluster_Amplitude *ReadEvent() { return p_ampl; }

    void SetAmpl(Cluster_Amplitude *ampl) { p_ampl = ampl; }

}; // end of class Achilles_Reader

} // end of namespace achilles
