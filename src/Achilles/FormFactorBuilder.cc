#include "Achilles/FormFactor.hh"
#include <iostream>

using achilles::FormFactorBuilder;

FormFactorBuilder &FormFactorBuilder::Vector(const std::string &name, const YAML::Node &node) {
    form_factor->vector = FormFactorFactory::Initialize(name, FFType::vector, node);
    return *this;
}

FormFactorBuilder &FormFactorBuilder::AxialVector(const std::string &name, const YAML::Node &node) {
    form_factor->axial = FormFactorFactory::Initialize(name, FFType::axial, node);
    return *this;
}

FormFactorBuilder &FormFactorBuilder::Coherent(const std::string &name, const YAML::Node &node) {
    form_factor->coherent = FormFactorFactory::Initialize(name, FFType::coherent, node);
    return *this;
}
