#pragma once

#include "METOOLS/Explicit/Form_Factor.H"
#include "METOOLS/Explicit/Vertex_Key.H"

#include "Achilles/FormFactor.hh"
#include "Achilles/Units.hh"

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
