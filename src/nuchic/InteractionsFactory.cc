#include "nuchic/Interactions.hh"

using namespace nuchic;

InteractionFactory& InteractionFactory::Instance() {
    static InteractionFactory factory = InteractionFactory();
    return factory;
}

bool InteractionFactory::Register(const std::string& name,
                                  InteractionFactory::TCreateMethod funcCreate) {
    auto it = methods.find(name);
    if(it == methods.end()) {
        methods.insert(std::make_pair(name, funcCreate));
        return true;
    }
    return false;
}

std::shared_ptr<Interactions> InteractionFactory::Create(const std::string& name,
                                                                         const std::string& filename) {
    auto it = methods.find(name);
    if(it != methods.end())
        return it -> second(filename);

    return nullptr;
}
