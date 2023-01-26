#include "Achilles/fortran/FNuclearModel.hh"
#include "Achilles/FourVector.hh"
#include "Achilles/Particle.hh"

using achilles::NuclearModel;
using achilles::FortranModel;

FortranModel::FortranModel(const YAML::Node &config, const YAML::Node &form_factor,
                           const std::shared_ptr<Nucleus> &nuc,
                           FormFactorBuilder &builder = FormFactorBuilder::Instance()) 
        : NuclearModel(form_factor, nuc, builder) {
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
    std::vector<NuclearModel::Currents> result(m_group.Processes().size());

    // Create momentum variables to pass to fortran
    const size_t nin = m_group.nucleons_in;
    const size_t nout = m_group.nucleons_out;
    std::vector<FourVector> moms;
    FourVector q;

    // Fill momentum
    // TODO: Make this cleaner, maybe pass FourVectors to fortran and use interface?
    // TODO: Make less dependent on 0th process
    int lepton = 1;
    auto ids = m_group.Processes()[0].Ids();
    for(size_t i = 0; i < ids.size(); ++i) {
        if(ParticleInfo(ids[i]).IsHadron()) {
            moms.push_back(event.Momentum()[i]);
        } else {
            // First lepton is always initial state, rest are final state
            q += lepton*event.Momentum()[i];
            lepton = -1;
        }
    }

    // Loop over the allowed initial states
    auto ffVals = EvalFormFactor(-q.M2()/1.0_GeV/1.0_GeV);
    for(size_t i = 0; i < m_group.Processes().size(); ++i) {
        auto process = m_group.Process(i);
        // Loop over the pid of q
        for(const auto &ffinfo : ff[i]) {
            // Load the form factors to pass to fortran
            // TODO: Have the model specify the needed form factors to load.
            //       This will be needed outside the interface as well
            auto formfactors = CouplingsFF(ffVals, ffinfo.second);

            // Get current from fortran and convert to right format
            std::complex<double> *cur = new std::complex<double>[NSpins()*4];
            size_t nff = formfactors.size();
            std::vector<long> pids_in;
            std::vector<long> pids_out;
            for(const auto &pid : process.state.first)
                pids_in.push_back(pid);
            for(const auto &pid : process.state.second)
                pids_out.push_back(pid);
            GetCurrents(pids_in.data(), pids_out.data(), moms.data(), nin, nout, &q,
                        formfactors.data(), nff, cur, NSpins(), 4);

            // Convert from array to Current
            for(size_t j = 0; j < NSpins(); ++j) {
                std::vector<std::complex<double>> tmp(4);
                for(size_t k = 0; k < tmp.size(); ++k) {
                    tmp[k] = cur[j + NSpins()*k];
                }
                result[i][ffinfo.first].push_back(tmp);
            }

            delete[] cur;
            cur = nullptr;
        }
    }

    return result;
}

std::unique_ptr<NuclearModel> FortranModel::Construct(const YAML::Node &config, const std::shared_ptr<Nucleus> &nuc) {
    auto form_factor = LoadFormFactor(config);
    return std::make_unique<FortranModel>(config, form_factor, nuc);
}
