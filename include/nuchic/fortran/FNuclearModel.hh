#ifndef FNUCLEAR_MODEL_HH
#define FNUCLEAR_MODEL_HH

#include "nuchic/NuclearModel.hh"

extern "C" {
    void RegisterAll();
    void ListModels();
    bool CreateModel(char*);
    bool InitModel(char*);
    void CleanUpEvent(std::complex<double>**, int*);
    void CleanUpModel();
    int GetMode();
    char* GetName();
    void GetCurrents(nuchic::FourVector *pin, size_t nin, size_t nout, nuchic::FourVector *q,
                     std::complex<double>* ff, size_t nff,
                     std::complex<double>** current, int *ncurrent);
    bool GetAllowedStates(nuchic::Process_Info*);
    size_t GetNSpins();
    bool FillNucleus(nuchic::Event*, const double*, size_t);
}

namespace nuchic {

class FortranModel : public NuclearModel, RegistrableNuclearModel<FortranModel> {
    public:
        FortranModel(const YAML::Node&, const YAML::Node&, FormFactorBuilder&);
        ~FortranModel() override { CleanUpModel(); }
        
        NuclearMode Mode() const override { 
            return static_cast<NuclearMode>(GetMode());
        }
        std::string PhaseSpace() const override {
            char *name = GetName();
            auto tmp = std::string(name);
            delete name;
            return tmp;
        }
        std::vector<Currents> CalcCurrents(const Event&, const std::vector<FFInfoMap>&) const override;
        void AllowedStates(Process_Info &info) override {
            GetAllowedStates(&info);
            m_info = info;
        }
        size_t NSpins() const override { return GetNSpins(); }
        bool FillNucleus(Event &event, const std::vector<double> &xsecs) const override { 
            return ::FillNucleus(&event, xsecs.data(), xsecs.size());
        }

        // Required factory methods
        static std::unique_ptr<NuclearModel> Construct(const YAML::Node&);
        static std::string Name() { return "FortranModel"; }

        // Method needed to register fortran models at start-up
        static void RegisterModels() { RegisterAll(); }
        static void DisplayModels() { 
            fmt::print("Available Fortran nuclear models:\n");
            ListModels();
        }

    private:
        Process_Info m_info;
};

}

#endif
