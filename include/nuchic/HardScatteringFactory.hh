#ifndef HARDSCATTERING_FACTORY_HH
#define HARDSCATTERING_FACTORY_HH

#include "nuchic/RunModes.hh"

#include <string> 
#include <memory>
#include <unordered_map>

namespace YAML {
class Node;
}

namespace nuchic {

class HardScattering;
class Beam;
class Nucleus;

class HardScatteringFactory {
    public:
        using TCreateMethod = std::unique_ptr<HardScattering>(*)(const YAML::Node&, RunMode);

        HardScatteringFactory() = delete;

        static bool Register(const std::string&, TCreateMethod);
        static std::unique_ptr<HardScattering> Create(const std::string&,
                const YAML::Node&, RunMode);
        static void ListHardScattering();

    private:
        static auto &methods() {
            static std::unordered_map<std::string, TCreateMethod> map;
            return map;
        }
};

#define REGISTER_HARDSCATTERING(hardscattering) \
    bool hardscattering::registered = nuchic::HardScatteringFactory::Register(hardscattering::GetName(), \
                                                                      hardscattering::Create)

}

#endif
