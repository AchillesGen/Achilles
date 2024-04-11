#ifndef FNUCLEAR_MODEL_HH
#define FNUCLEAR_MODEL_HH

#include "Achilles/NuclearModel.hh"

extern "C" {
void RegisterAll();
void ListModels();
bool CreateModel(char *);
bool InitModel(char *);
void CleanUpEvent(std::complex<double> **, int *);
void CleanUpModel();
int GetMode();
int GetFrame();
char *ModelName();
char *GetName_();
void GetCurrents(long *pids_in, long *pids_out, achilles::FourVector *pin, size_t nin, size_t nout,
                 const achilles::FourVector *q, std::map<std::string, std::complex<double>> *ff,
                 std::complex<double> *current, size_t nspins, size_t nlorentz);
double GetInitialStateWeight(long *pids_in, achilles::FourVector *pin, size_t nin, size_t nproton,
                             size_t nneutron);
}

namespace achilles {

class FortranModel : public NuclearModel, RegistrableNuclearModel<FortranModel> {
  public:
    FortranModel(const YAML::Node &, const YAML::Node &, FormFactorBuilder &);
    ~FortranModel() override { CleanUpModel(); }

    NuclearMode Mode() const override { return static_cast<NuclearMode>(GetMode()); }
    std::string PhaseSpace(PID) const override;
    NuclearFrame Frame() const override { 
        return static_cast<NuclearFrame>(GetFrame());
    }
    Currents CalcCurrents(const std::vector<Particle> &, const std::vector<Particle> &,
                          const FourVector &, const FFInfoMap &) const override;
    double InitialStateWeight(const std::vector<Particle> &, size_t, size_t) const override;

    // Required factory methods
    static std::unique_ptr<NuclearModel> Construct(const YAML::Node &);
    static std::string Name() { return "FortranModel"; }

    // TODO: Allow fortran codes to fill these out
    std::string GetName() const override { return ModelName(); }
    std::string InspireHEP() const override { return ""; }

    // Method needed to register fortran models at start-up
    static void RegisterModels() { RegisterAll(); }
    static void DisplayModels() {
        fmt::print("Available Fortran nuclear models:\n");
        ListModels();
    }

  private:
    mutable bool is_hydrogen;
    const WardGauge m_ward;
};

} // namespace achilles

#endif
