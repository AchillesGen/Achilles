#include "Achilles/fortran/FNuclearModel.hh"
#include "Achilles/FourVector.hh"
#include "Achilles/Particle.hh"
#include "Achilles/Poincare.hh"

using achilles::FortranModel;
using achilles::NuclearModel;

FortranModel::FortranModel(const YAML::Node &config, const YAML::Node &form_factor,
                           FormFactorBuilder &builder = FormFactorBuilder::Instance())
    : NuclearModel(form_factor, builder),
    m_ward{ToEnum(config["NuclearModel"]["Ward"].as<std::string>())}
    {
    // Setup fortran side
    auto modelname = config["NuclearModel"]["Name"].as<std::string>();
    size_t len_model = modelname.size();
    auto cmodelname = std::make_unique<char[]>(len_model + 1);
    strcpy(cmodelname.get(), modelname.c_str());

    if(!CreateModel(cmodelname.get())) {
        auto msg = fmt::format(
            "NuclearModel: Invalid model requested {}. Please ensure the model is registered",
            modelname);
        throw std::runtime_error(msg);
    }

    auto filename = config["NuclearModel"]["ConfigFile"].as<std::string>();
    size_t len = filename.size();
    auto cfilename = std::make_unique<char[]>(len + 1);
    strcpy(cfilename.get(), filename.c_str());

    if(!InitModel(cfilename.get())) {
        auto msg = fmt::format("NuclearModel: Could not initialize model {} using file {}.",
                               modelname, filename);
        throw std::runtime_error(msg);
    }
}

NuclearModel::Currents FortranModel::CalcCurrents(const std::vector<Particle> &had_in,
                                                  const std::vector<Particle> &had_out,
                                                  const FourVector &qVec,
                                                  const FFInfoMap &ff) const {

    
    if(had_in[0].ID() == PID::neutron() && is_hydrogen) return {};
    
    NuclearModel::Currents result;

    // Create momentum variables to pass to fortran
    const size_t nin = had_in.size();
    const size_t nout = had_out.size();
    spdlog::debug("q = {}", qVec);

    auto omega = qVec[0];

    auto ffVals = EvalFormFactor(-qVec.M2() / 1.0_GeV / 1.0_GeV);
    std::vector<long> pids_in;
    std::vector<long> pids_out;
    std::vector<FourVector> moms;
    for(const auto &part : had_in) {
        pids_in.push_back(part.ID());
        moms.push_back(part.Momentum());
    }
    for(const auto &part : had_out) {
        pids_out.push_back(part.ID());
        moms.push_back(part.Momentum());
    }

    // Loop over the pid of q
    for(const auto &ffinfo : ff) {
        // Load the form factors to pass to fortran
        // TODO: Have the model specify the needed form factors to load.
        //       This will be needed outside the interface as well
        auto formfactors = CouplingsFF(ffVals, ffinfo.second);

        // Get current from fortran and convert to right format
        std::complex<double> *cur = new std::complex<double>[NSpins() * 4];

        std::map<std::string, std::complex<double>> ffmap;
        for(const auto &factor : formfactors) { ffmap[ToString(factor.first)] = factor.second; }

        GetCurrents(pids_in.data(), pids_out.data(), moms.data(), nin, nout, &qVec, &ffmap, cur,
                    NSpins(), 4);

        // Convert from array to Current
        for(size_t j = 0; j < NSpins(); ++j) {
            std::vector<std::complex<double>> tmp(4);
            for(size_t k = 0; k < tmp.size(); ++k) { tmp[k] = cur[j + NSpins() * k]; }
            result[ffinfo.first].push_back(tmp);
        }

        // Correct the Ward identity
        for(auto &subcur : result[ffinfo.first]) {
            switch(m_ward) {
            case WardGauge::None:
                continue;
                break;
            case WardGauge::Coulomb:
                CoulombGauge(subcur, qVec, omega);
                break;
            case WardGauge::Weyl:
                WeylGauge(subcur, qVec, omega);
                break;
            case WardGauge::Landau:
                LandauGauge(subcur, qVec);
                break;
            }
        }

        delete[] cur;
        cur = nullptr;
    }

    return result;
}

void FortranModel::CoulombGauge(VCurrent &cur, const FourVector &q, double omega) const {
    FourVector cur4_real{cur[0].real(), cur[1].real(), cur[2].real(), cur[3].real()};
    FourVector cur4_imag{cur[0].imag(), cur[1].imag(), cur[2].imag(), cur[3].imag()};
    FourVector ref{0, 0, 0, 1};
    Poincare poincare(q, ref, 0);
    poincare.Rotate(cur4_real);
    poincare.Rotate(cur4_imag);

    cur4_real[3] = omega / q.P() * cur[0].real();
    cur4_imag[3] = omega / q.P() * cur[0].imag();
    poincare.RotateBack(cur4_real);
    poincare.RotateBack(cur4_imag);
    for(size_t i = 0; i < 4; ++i) cur[i] = {cur4_real[i], cur4_imag[i]};
}

void FortranModel::WeylGauge(VCurrent &cur, const FourVector &q, double omega) const {
    FourVector cur4_real{cur[0].real(), cur[1].real(), cur[2].real(), cur[3].real()};
    FourVector cur4_imag{cur[0].imag(), cur[1].imag(), cur[2].imag(), cur[3].imag()};
    FourVector ref{0, 0, 0, 1};
    Poincare poincare(q, ref, 0);

    poincare.Rotate(cur4_real);
    poincare.Rotate(cur4_imag);

    cur4_real[0] = q.P() / omega * cur4_real[3];
    cur4_imag[0] = q.P() / omega * cur4_imag[3];
    poincare.RotateBack(cur4_real);
    poincare.RotateBack(cur4_imag);
    for(size_t i = 0; i < 4; ++i) cur[i] = {cur4_real[i], cur4_imag[i]};
}

void FortranModel::LandauGauge(VCurrent &cur, const FourVector &q) const {
    auto jdotq = cur[0] * q[0] - cur[1] * q[1] - cur[2] * q[2] - cur[3] * q[3];
    auto Q2 = -q.M2();
    for(size_t i = 0; i < 4; ++i) cur[i] += jdotq / Q2 * q[i];
}

double FortranModel::InitialStateWeight(const std::vector<Particle> &had_in, size_t nproton,
                                        size_t nneutron) const {

    if(is_hydrogen) return had_in[0].ID() == PID::proton() ? 1 : 0;

    size_t nin = had_in.size();
    std::vector<long> pids_in;
    std::vector<FourVector> moms;
    for(const auto &part : had_in) {
        pids_in.push_back(part.ID());
        moms.push_back(part.Momentum());
    }

    return GetInitialStateWeight(pids_in.data(), moms.data(), nin, nproton, nneutron);
}

std::unique_ptr<NuclearModel> FortranModel::Construct(const YAML::Node &config) {
    auto form_factor = LoadFormFactor(config);
    return std::make_unique<FortranModel>(config, form_factor);
}

std::string FortranModel::PhaseSpace(PID nuc_pid) const {
    if(nuc_pid != PID::hydrogen()){
        char *name = GetName_();
        auto tmp = std::string(name);
        delete name;
        return tmp;
    }
    is_hydrogen = true;
    return "Coherent";
}
