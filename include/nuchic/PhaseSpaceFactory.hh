#ifndef PHASE_SPACE_FACTORY_HH
#define PHASE_SPACE_FACTORY_HH

#include "spdlog/spdlog.h"

#include <functional>
#include <map>
#include <memory>

namespace nuchic {

template<class Base, typename T>
class PSFactory {
    using Constructor = std::function<std::unique_ptr<Base>(T)>;

    static std::map<std::string, Constructor>& Registry() {
        static std::map<std::string, Constructor> registry;
        return registry;
    }

    public:
        static std::unique_ptr<Base> Build(const std::string &name, const T &t) {
            auto constructor = Registry().at(name);
            return constructor(t);
        }

        template<class Derived>
        static void Register(std::string name) {
            if(IsRegistered(name))
                spdlog::error("{} is already registered!", name);
            spdlog::trace("Registering {} Phase Space", name);
            Registry()[name] = Derived::Construct;
        }

        static bool IsRegistered(std::string name) {
            return Registry().find(name) != Registry().end();
        }

        static void DisplayPhaseSpaces() {
            fmt::print("Registered {} Phase Spaces:\n", Base::Name());
            for(const auto &registered : Registry())
                fmt::print("  - {}\n", registered.first);
        }
};

template<class Base, class Derived, typename T>
class RegistrablePS {
    protected:
        RegistrablePS() = default;
        virtual ~RegistrablePS() {
            if(!m_registered)
                spdlog::error("Error registering phase space");
        }

        static bool Register() {
            PSFactory<Base, T>::template Register<Derived>(Derived::Name());
            return true;
        }

    private:
        static const bool m_registered;
};
template<class Base, class Derived, typename T>
const bool RegistrablePS<Base, Derived, T>::m_registered = RegistrablePS<Base, Derived, T>::Register();

}

#endif
