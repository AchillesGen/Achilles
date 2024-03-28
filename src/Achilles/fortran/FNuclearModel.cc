#include "Achilles/fortran/FNuclearModel.hh"
#include "Achilles/FourVector.hh"
#include "Achilles/Particle.hh"

using achilles::FortranModel;
using achilles::NuclearModel;

FortranModel::FortranModel(const YAML::Node &config, const YAML::Node &form_factor,
                           FormFactorBuilder &builder = FormFactorBuilder::Instance())
    : NuclearModel(form_factor, builder) {
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
    NuclearModel::Currents result;

    // Create momentum variables to pass to fortran
    const size_t nin = had_in.size();
    const size_t nout = had_out.size();
    spdlog::debug("q = {}", qVec);

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

        delete[] cur;
        cur = nullptr;
    }

    return result;
}

double FortranModel::InitialStateWeight(const std::vector<Particle> &had_in, size_t nproton,
                                        size_t nneutron) const {
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
