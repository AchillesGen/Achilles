#include "plugins/Sherpa/FormFactors.hh"

#include "ATOOLS/Org/Message.H"
#include "Achilles/Constants.hh"
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
    // TODO: Fix sign in Model file
    // TODO: Include factor of 1/(2Mn) for F2 in Model file
    FormFactor::Values ffs((*p_ff)(-Q2));
    if(m_id == "F1p")
        result = ffs.F1p;
    else if(m_id == "F2p")
        result = ffs.F2p / (2 * Constant::mN / 1_GeV);
    else if(m_id == "F1n")
        result = ffs.F1n;
    else if(m_id == "F2n")
        result = ffs.F2n / (2 * Constant::mN / 1_GeV);
    else if(m_id == "FA")
        result = ffs.FA;
    else if(m_id == "FCoh")
        result = ffs.Fcoh;
    else
        THROW(not_implemented, "Form factor not found");
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

FORM_FACTOR_IMPLEMENTATION(F1p)
FORM_FACTOR_IMPLEMENTATION(F2p)
FORM_FACTOR_IMPLEMENTATION(F1n)
FORM_FACTOR_IMPLEMENTATION(F2n)
FORM_FACTOR_IMPLEMENTATION(FA)
FORM_FACTOR_IMPLEMENTATION(FCoh)

#undef FORM_FACTOR_IMPLEMENTATION
