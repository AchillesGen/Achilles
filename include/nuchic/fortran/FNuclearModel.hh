#ifndef FNUCLEAR_MODEL_HH
#define FNUCLEAR_MODEL_HH

#include "nuchic/NuclearModel.hh"

extern "C" {
    void RegisterAll();
    bool CreateModel(char*);
    bool InitModel(char*);
    void CleanUpEvent(std::complex<double>**, size_t*);
    void CleanUpModel();
    int GetMode();
    void GetName(char*);
    // TODO: Figure out how to pass a FFInfoMap to fortran
    //       - Evaluate form factors in C++ and just pass the complex valued array to fortran
    //       - Develop an interface to the form factor and maps
    //       - Other ideas?
    void GetCurrents(const nuchic::Event*, int**, size_t*, std::complex<double>**, size_t*);
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
            char *name = nullptr;
            GetName(name);
            auto tmp = std::string(name);
            delete name;
            return tmp;
        }
        std::vector<Currents> CalcCurrents(const Event&, const std::vector<FFInfoMap>&) const override;
        void AllowedStates(Process_Info &info) const override {
            GetAllowedStates(&info);
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
};

}

#endif
