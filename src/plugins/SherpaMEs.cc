#include "plugins/SherpaMEs.hh"
#include "SHERPA/Main/Sherpa.H"
#include "SHERPA/Initialization/Initialization_Handler.H"
#include "SHERPA/PerturbativePhysics/Matrix_Element_Handler.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "COMIX/Main/Single_Process.H"
#include "AMEGIC++/Main/Process_Base.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "METOOLS/Explicit/Color_Calculator.H"
#include "ATOOLS/Org/Library_Loader.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/My_File.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "nuchic/Logging.hh"

using namespace SHERPA;
using namespace PHASIC;
using namespace ATOOLS;
using namespace nuchic;

SherpaMEs::~SherpaMEs()
{
  if (p_sherpa) delete p_sherpa;
}

void SherpaMEs::addParameter(std::vector<char*>& argv,const std::string& val) const
{
  spdlog::info("ShowerMEsSherpa::init: Setting parameter '{}'",val);
  argv.push_back(new char[val.length()+1]);
  strcpy(argv.back(),val.c_str());
}

int SherpaMEs::SherpaVerbosity(int loglevel) const
{
  switch(loglevel) {
  case 0/*trace*/: return 47;
  case 1/*debug*/: return 15;
  case 2/*info*/: return 2;
  case 3/*warn*/: return 1;
  case 4/*error*/: return 1;
  case 5/*critical*/: return 0;
  }
  return 0;
}

bool SherpaMEs::Initialize(const std::vector<std::string> &args)
{
  p_sherpa = new Sherpa();
  std::vector<char*> argv;
  addParameter(argv,"Sherpa");
  addParameter(argv,"EVENTS=0");
  addParameter(argv,"INIT_ONLY=6");
  int level(SherpaVerbosity(spdlog::get("nuchic")->level()));
  addParameter(argv,"OUTPUT="+ToString(level));
  // Set beams and PDFs.
  addParameter(argv,"BEAM_1=90");
  addParameter(argv,"BEAM_2=90");
  addParameter(argv,"BEAM_ENERGY_1=100");
  addParameter(argv,"BEAM_ENERGY_2=100");
  addParameter(argv,"PDF_SET=None");
  addParameter(argv,"PDF_LIBRARY=None");
  addParameter(argv,"ME_SIGNAL_GENERATOR=Comix");
  // add additional commandline parameters
  for (const auto &arg: args) addParameter(argv,arg);
  // Initialise Sherpa and return.
  p_sherpa->InitializeTheRun(argv.size(),&argv[0]);
  return true;
}

Process_Base* SherpaMEs::getProcess(Cluster_Amplitude* const ampl) {
  // Initialise process map.
  std::string megen = "Comix";
  // StringProcess_Map* pm(m_pmap[nlo_type::lo]);
  My_In_File::OpenDB(rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Sherpa/");
  My_In_File::OpenDB(rpa->gen.Variable("SHERPA_CPP_PATH")
    +"/Process/"+megen+"/");

  // Initialise process info and create label.
  PHASIC::Process_Info pi;
  pi.m_maxcpl[0]   = pi.m_mincpl[0] = ampl->OrderQCD();
  pi.m_maxcpl[1]   = pi.m_mincpl[1] = ampl->OrderEW();
  pi.m_megenerator = megen;
  for (size_t i(0); i<ampl->NIn(); ++i) {
    Flavour fl(ampl->Leg(i)->Flav().Bar());
    if (Flavour(kf_jet).Includes(fl)) fl=Flavour(kf_jet);
    pi.m_ii.m_ps.push_back(Subprocess_Info(fl,"",""));
  }
  for (size_t i(ampl->NIn()); i<ampl->Legs().size(); ++i) {
    Flavour fl(ampl->Leg(i)->Flav());
    if (Flavour(kf_jet).Includes(fl)) fl=Flavour(kf_jet);
    pi.m_fi.m_ps.push_back(Subprocess_Info(fl,"",""));
  }

  // Create and initialise process.
  Process_Base* proc = p_sherpa->GetInitHandler()->
    GetMatrixElementHandler()->Generators()->InitializeProcess(pi,false);
  if (proc == nullptr) {
    My_In_File::CloseDB(rpa->gen.Variable("SHERPA_CPP_PATH")
      +"/Process/"+megen+"/");
    My_In_File::CloseDB(rpa->gen.Variable("SHERPA_CPP_PATH")
      +"/Process/Sherpa/");
    return nullptr;
  }
  // m_procs.push_back(proc);
  proc->SetSelector(Selector_Key(nullptr,new Data_Reader(),true));
  proc->SetScale(Scale_Setter_Arguments(MODEL::s_model,
    "VAR{100}{100}", "Alpha_QCD 1"));
  proc->SetKFactor(KFactor_Setter_Arguments("NO"));
  if (proc->Get<COMIX::Process_Base>()) {
    proc->Get<COMIX::Process_Base>()->Tests();
  }
  // proc->FillProcessMap(&m_pmap);
  proc->ConstructColorMatrix();
  My_In_File::CloseDB(rpa->gen.Variable("SHERPA_CPP_PATH")
    +"/Process/"+megen+"/");
  My_In_File::CloseDB(rpa->gen.Variable("SHERPA_CPP_PATH")
    +"/Process/Sherpa/");

  return proc;
}

bool SherpaMEs::InitializeProcess(const Process_Info &info)
{
  Cluster_Amplitude* ampl = Cluster_Amplitude::New();
  int nqcd(0), nIn(0), cmin(std::numeric_limits<int>::max()), cmax(0);
  for (size_t i(0);i<info.m_ids.size();++i) {
    ampl->CreateLeg(Vec4D(),i<2?Flavour(info.m_ids[i]):Flavour(info.m_ids[i]).Bar());
  }
  ampl->SetNIn(2);
  ampl->SetOrderQCD(0);
  ampl->SetOrderEW(info.m_ids.size()-2);
  Process_Base::SortFlavours(ampl);
  // StringProcess_Map* pm = m_pmap[nlo_type::lo];
  std::string name = Process_Base::GenerateName(ampl);
  // if (pm->find(name) == pm->end()) {
  //   // infoPtr->errorMsg("Could not find process "+name);
  //   getProcess(ampl);
  // }
  // Process_Base* proc = pm->find(name)->second;
  // if (proc == nullptr) {
  //   // infoPtr->errorMsg("Error in ShowerMEsSherpa: process not found.");
  //   return false;
  // }
  return true;
}
