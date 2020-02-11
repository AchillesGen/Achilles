#include "nuchic/Interactions.hh"

nuchic::InteractionFactory& nuchic::InteractionFactory::Instance() {
    static nuchic::InteractionFactory factory = InteractionFactory();
    return factory;
}

bool nuchic::InteractionFactory::Register(const std::string& name,
                                  nuchic::InteractionFactory::TCreateMethod funcCreate) {
    auto it = methods.find(name);
    if(it == methods.end()) {
        methods.insert(std::make_pair(name, funcCreate));
        return true;
    }
    return false;
}

std::shared_ptr<nuchic::Interactions> nuchic::InteractionFactory::Create(const std::string& name,
                                                                         const std::string& filename) {
    auto it = methods.find(name);
    if(it != methods.end())
        return it -> second(filename);

    return nullptr;
}
