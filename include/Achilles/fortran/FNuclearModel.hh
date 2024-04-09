#ifndef FNUCLEAR_MODEL_HH
#define FNUCLEAR_MODEL_HH

#include "Achilles/NuclearModel.hh"

extern "C" {
void RegisterAll();
void ListModels();
bool CreateModel(char *, size_t &);
bool InitModel(char *, std::map<std::string, double> *pm, size_t);
void CleanUpEvent(std::complex<double> **, int *);
void CleanUpModel(size_t);
int GetMode(size_t);
int GetFrame(size_t);
char *ModelName(size_t);
char *ModelPS(size_t);
void GetCurrents(size_t, long *pids_in, long *pids_out, long *pids_spect, achilles::FourVector *pin,
                 size_t nin, size_t nout, size_t nspect, const achilles::FourVector *q,
                 std::map<std::string, std::complex<double>> *ff, std::complex<double> *current,
                 size_t nspins, size_t nlorentz);
double GetInitialStateWeight(size_t, long *pids_in, long *pids_spect, achilles::FourVector *mom,
                             size_t nin, size_t nspect, size_t nproton, size_t nneutron);
char *GetInspireHEP(size_t);
}

namespace achilles {

class FortranModel : public NuclearModel, RegistrableNuclearModel<FortranModel> {
  public:
    FortranModel(const YAML::Node &, const YAML::Node &, FormFactorBuilder &);
    ~FortranModel() override { CleanUpModel(m_model); }

    NuclearMode Mode() const override { return static_cast<NuclearMode>(GetMode(m_model)); }
    std::string PhaseSpace(PID) const override;
    NuclearFrame Frame() const override { return static_cast<NuclearFrame>(GetFrame(m_model)); }
    Currents CalcCurrents(const std::vector<Particle> &, const std::vector<Particle> &,
                          const std::vector<Particle> &, const FourVector &,
                          const FFInfoMap &) const override;
    double InitialStateWeight(const std::vector<Particle> &, const std::vector<Particle> &, size_t,
                              size_t) const override;

    // Required factory methods
    static std::unique_ptr<NuclearModel> Construct(const YAML::Node &);
    static std::string Name() { return "FortranModel"; }
    std::string PSName() const override { return ModelPS(m_model); }

    std::string GetName() const override { return ModelName(m_model); }
    std::string InspireHEP() const override { return GetInspireHEP(m_model); }

    // Method needed to register fortran models at start-up
    static void RegisterModels() { RegisterAll(); }
    static void DisplayModels() {
        fmt::print("Available Fortran nuclear models:\n");
        ListModels();
    }

  private:
    mutable bool is_hydrogen{false};
    mutable bool is_free_neutron{false};
    const WardGauge m_ward;
    size_t m_model;
    std::map<std::string, double> param_map;
};

} // namespace achilles

#endif
