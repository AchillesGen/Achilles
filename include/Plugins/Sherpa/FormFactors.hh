#pragma once

#include "Achilles/FormFactor.hh"

#ifdef UFO_2_0
#include "METOOLS/Explicit/Form_Factor.H"
#include "METOOLS/Explicit/Vertex_Key.H"

namespace achilles {

class FormFactorInterface : public METOOLS::Form_Factor {
  private:
    static std::unique_ptr<FormFactor> p_ff;
    std::string m_id;
    int m_mode;

  public:
    FormFactorInterface(const METOOLS::Vertex_Key &key, const std::string &id);
    double FF();
    static void SetFormFactor(std::unique_ptr<FormFactor> ff) { p_ff = std::move(ff); }
}; // end of class FormFactorInterface

} // end of namespace achilles

#else

namespace achilles {

class FormFactorInterface {
  public:
    static void SetFormFactor(std::unique_ptr<FormFactor>) {
        throw std::runtime_error("FormFactorInterface::SetFormFactor not implemented");
    }
};
} // end of namespace achilles

#endif
