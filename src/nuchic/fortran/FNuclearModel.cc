#include "nuchic/fortran/FNuclearModel.hh"

using nuchic::NuclearModel;
using nuchic::FortranModel;

FortranModel::FortranModel(const YAML::Node &config, const YAML::Node &form_factor, FormFactorBuilder &builder = FormFactorBuilder::Instance()) 
        : NuclearModel(form_factor, builder) {
    // Setup fortran side
    auto modelname = config["Name"].as<std::string>();
    size_t len_model = modelname.size();
    auto cmodelname = std::make_unique<char[]>(len_model+1);
    strcpy(cmodelname.get(), modelname.c_str());

    if(!CreateModel(cmodelname.get())) {
        auto msg = fmt::format("NuclearModel: Invalid model requested {}. Please ensure the model is registered", modelname);
        throw std::runtime_error(msg);
    }

    auto filename = config["configFile"].as<std::string>();
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
    std::vector<NuclearModel::Currents> result;
    std::complex<double> *cur = nullptr;
    std::vector<int> bosons;
    int nbosons{};
    size_t size = 4*GetNSpins();

    // Get current from fortran and convert to right format
    // TODO: How to handle various initial states?
    GetCurrents(&event, bosons.data(), &nbosons, &cur, &size);
    // TODO: Figure out how to handle converting to std::vector<Currents>

    // Clean up memory usage
    CleanUpEvent(&cur, &size);
    cur = nullptr;

    return result;
}

std::unique_ptr<NuclearModel> FortranModel::Construct(const YAML::Node &config) {
    auto form_factor = LoadFormFactor(config);
    return std::make_unique<FortranModel>(config, form_factor);
}
