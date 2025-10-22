#include "Plugins/Sherpa/FormFactors.hh"
#include "Achilles/Constants.hh"

#include "ATOOLS/Org/Message.H"
#include "Achilles/Constants.hh"
#include "Achilles/ParticleInfo.hh"
#include "METOOLS/Explicit/Current.H"
#include "METOOLS/Explicit/Vertex.H"
#include "MODEL/Main/Model_Base.H"

#include "ATOOLS/Org/Exception.H"

using namespace ATOOLS;
using namespace METOOLS;
using namespace achilles;

std::shared_ptr<FormFactor> FormFactorInterface::p_ff(nullptr);

FormFactorInterface::FormFactorInterface(const Vertex_Key &key, const std::string &id)
    : Form_Factor(id, key), m_id(id), m_mode(-1) {
    for(size_t i(0); i < key.m_j.size(); ++i)
        if(key.m_j[i]->Flav().IsPhoton()) m_mode = static_cast<int>(i);
}

double FormFactorInterface::FF(double Q2) const {
    double result(0.);
    FormFactor::Values ffs((*p_ff)(Q2));
    if(m_id == "F1PA" || m_id == "F1PZ" || m_id == "F1PW") {
        result = ffs.F1p;
        // spdlog::trace("m_id = {}, Q2 = {}, F1p = {}", m_id, Q2, result);
    } else if(m_id == "F2PA" || m_id == "F2PZ" || m_id == "F2PW") {
        result = ffs.F2p;
        // spdlog::trace("m_id = {}, Q2 = {}, F2p = {}", m_id, Q2, result);
    } else if(m_id == "F1NA" || m_id == "F1NZ" || m_id == "F1NW") {
        result = ffs.F1n;
        // spdlog::trace("m_id = {}, Q2 = {}, F1n = {}", m_id, Q2, result);
    } else if(m_id == "F2NA" || m_id == "F2NZ" || m_id == "F2NW") {
        result = ffs.F2n;
        // spdlog::trace("m_id = {}, Q2 = {}, F2n = {}", m_id, Q2, result);
    } else if(m_id == "FA") {
        result = ffs.FA;
        // spdlog::trace("m_id = {}, Q2 = {}, FA = {}", m_id, Q2, result);
    } else if(m_id == "FCoh")
        result = ffs.Fcoh;
    else if(m_id == "FPiA") {
        result = ffs.FPiA;
        // spdlog::trace("Q2 = {}, FPiA = {}", Q2, result);
    } else if(m_id == "FPiW")
        result = ffs.FPiW;
    else if(m_id == "FPiZ")
        result = ffs.FPiZ;
    else if(m_id == "FGAPi") {
        result = ffs.FGAPi;
        // spdlog::trace("Q2 = {}, FGAPi = {}", Q2, result);
    } else if(m_id == "FDeltaPG1") {
        result = ffs.Gep + ffs.Gmp;
    } else if(m_id == "FDeltaPG2") {
        double MNucleon = achilles::ParticleInfo(achilles::PID::proton()).Mass();
        double MDelta = achilles::ParticleInfo(achilles::PID::deltap()).Mass();
        result = 2 * ffs.Gep * (MNucleon + MDelta) / (MNucleon - MDelta);
    } else if(m_id == "FDeltaNG1") {
        result = ffs.Gen + ffs.Gmn;
    } else if(m_id == "FDeltaNG2") {
        double MNucleon = achilles::ParticleInfo(achilles::PID::neutron()).Mass();
        double MDelta = achilles::ParticleInfo(achilles::PID::delta0()).Mass();
        result = 2 * ffs.Gen * (MNucleon + MDelta) / (MNucleon - MDelta);
    } else
        THROW(not_implemented, "Form factor not found");

    spdlog::trace("{} = {} (Gep = {})", m_id, result, ffs.Gep);
    return result;
}

#define FORM_FACTOR_IMPLEMENTATION(NAME)                                                         \
    namespace achilles {                                                                         \
    class NAME : public Getter_Function<Form_Factor, Vertex_Key, std::less<std::string> > {      \
        static NAME s_initializer;                                                               \
                                                                                                 \
      protected:                                                                                 \
        void PrintInfo(std::ostream &str, const size_t) const {                                  \
            str << #NAME;                                                                        \
        }                                                                                        \
        Object_Type *operator()(const Parameter_Type &args) const {                              \
            return new FormFactorInterface(args, #NAME);                                         \
        }                                                                                        \
                                                                                                 \
      public:                                                                                    \
        NAME(const std::string &name)                                                            \
            : ATOOLS::Getter_Function<Form_Factor, Vertex_Key, std::less<std::string> >(name) {} \
    };                                                                                           \
    NAME NAME::s_initializer(#NAME);                                                             \
    }

FORM_FACTOR_IMPLEMENTATION(F1PA)
FORM_FACTOR_IMPLEMENTATION(F2PA)
FORM_FACTOR_IMPLEMENTATION(F1NA)
FORM_FACTOR_IMPLEMENTATION(F2NA)
FORM_FACTOR_IMPLEMENTATION(F1PZ)
FORM_FACTOR_IMPLEMENTATION(F2PZ)
FORM_FACTOR_IMPLEMENTATION(F1NZ)
FORM_FACTOR_IMPLEMENTATION(F2NZ)
FORM_FACTOR_IMPLEMENTATION(F1PW)
FORM_FACTOR_IMPLEMENTATION(F2PW)
FORM_FACTOR_IMPLEMENTATION(F1NW)
FORM_FACTOR_IMPLEMENTATION(F2NW)
FORM_FACTOR_IMPLEMENTATION(FA)
FORM_FACTOR_IMPLEMENTATION(FCoh)
FORM_FACTOR_IMPLEMENTATION(FPiA)
FORM_FACTOR_IMPLEMENTATION(FPiW)
FORM_FACTOR_IMPLEMENTATION(FPiZ)
FORM_FACTOR_IMPLEMENTATION(FGAPi)
FORM_FACTOR_IMPLEMENTATION(FDeltaPG1)
FORM_FACTOR_IMPLEMENTATION(FDeltaPG2)
FORM_FACTOR_IMPLEMENTATION(FDeltaNG1)
FORM_FACTOR_IMPLEMENTATION(FDeltaNG2)

#undef FORM_FACTOR_IMPLEMENTATION
