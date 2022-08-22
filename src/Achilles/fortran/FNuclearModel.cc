#include "Achilles/fortran/FNuclearModel.hh"
#include "Achilles/FourVector.hh"
#include "Achilles/Particle.hh"

using achilles::NuclearModel;
using achilles::FortranModel;

FortranModel::FortranModel(const YAML::Node &config, const YAML::Node &form_factor, FormFactorBuilder &builder = FormFactorBuilder::Instance()) 
        : NuclearModel(form_factor, builder) {
    // Setup fortran side
    auto modelname = config["NuclearModel"]["Name"].as<std::string>();
    size_t len_model = modelname.size();
    auto cmodelname = std::make_unique<char[]>(len_model+1);
    strcpy(cmodelname.get(), modelname.c_str());

    if(!CreateModel(cmodelname.get())) {
        auto msg = fmt::format("NuclearModel: Invalid model requested {}. Please ensure the model is registered", modelname);
        throw std::runtime_error(msg);
    }

    auto filename = config["NuclearModel"]["ConfigFile"].as<std::string>();
    size_t len = filename.size();
    auto cfilename = std::make_unique<char[]>(len+1);
    strcpy(cfilename.get(), filename.c_str());

    if(!InitModel(cfilename.get())) {
        auto msg = fmt::format("NuclearModel: Could not initialize model {} using file {}.",
                               modelname, filename);
        throw std::runtime_error(msg);
    }
}

std::vector<NuclearModel::Currents> FortranModel::CalcCurrents(const Event &event,
                                                               const std::vector<FFInfoMap> &ff) const {
    std::vector<NuclearModel::Currents> result(m_info.m_states.size());

    // Create momentum variables to pass to fortran
    const size_t nin = m_info.m_states.begin()->first.size();
    const size_t nout = m_info.m_states.begin()->second.size();
    std::vector<FourVector> moms;
    FourVector q;

    // Fill momentum
    // TODO: Make this cleaner, maybe pass FourVectors to fortran and use interface?
    int lepton = 1;
    auto ids = m_info.Ids();
    for(size_t i = 0; i < ids.size(); ++i) {
        if(ParticleInfo(ids[i]).IsHadron()) {
            moms.push_back(event.Momentum()[i]);
        } else {
            // First lepton is always initial state, rest are final state
            q += lepton*event.Momentum()[i];
            lepton = -1;
        }
    }

    // Ensure the momentum vector is the right size
    if(moms.size() != nin + nout)
        throw std::runtime_error("FortranModel: Invalid number of hadrons in event.");

    // Loop over the allowed initial states
    auto ffVals = EvalFormFactor(-q.M2()/1.0_GeV/1.0_GeV);
    for(size_t i = 0; i < m_info.m_states.size(); ++i) {
        // Loop over the pid of q
        for(const auto &ffinfo : ff[i]) {
            // Load the form factors to pass to fortran
            // TODO: Have the model specify the needed form factors to load.
            //       This will be needed outside the interface as well
            auto formfactors = CouplingsFF(ffVals, ffinfo.second);

            // Get current from fortran and convert to right format
            std::complex<double> *cur = nullptr;
            size_t nff = formfactors.size();
            int ncur{};
            GetCurrents(moms.data(), nin, nout, &q, formfactors.data(), nff, &cur, &ncur);

            // Convert from array to Current
            for(size_t j = 0; j < static_cast<size_t>(ncur/4); ++j) {
                std::vector<std::complex<double>> tmp(4);
                for(size_t mu = 0; mu < 4; ++mu) {
                    tmp[mu] = cur[mu+4*j];
                }
                result[i][ffinfo.first].push_back(tmp);
            }

            // Clean up memory usage
            CleanUpEvent(&cur, &ncur);
            cur = nullptr;
        }
    }

    return result;
}

std::unique_ptr<NuclearModel> FortranModel::Construct(const YAML::Node &config) {
    auto form_factor = LoadFormFactor(config);
    return std::make_unique<FortranModel>(config, form_factor);
}
