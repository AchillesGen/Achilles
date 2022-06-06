#ifndef FACTORY_HH
#define FACTORY_HH

#include "spdlog/spdlog.h"
#include <functional>
#include <map>
#include <memory>

namespace achilles {

template<typename Base, typename ...Args>
class Factory {
    using Constructor = std::function<std::unique_ptr<Base>(Args&&...)>;

    static std::map<std::string, Constructor>& Registry() {
        static std::map<std::string, Constructor> registry;
        return registry;
    }

    public:
        static std::unique_ptr<Base> Initialize(const std::string &name,
                                                Args &&...args) {
            auto constructor = Registry().at(name);
            return constructor(std::forward<Args>(args)...);
        }

        template<class Derived>
        static void Register(const std::string &name) {
            if(IsRegistered(name)) 
                spdlog::error("{} is already registered!", name);
            spdlog::trace("Registering {}", name);
            Registry()[name] = Derived::Construct;
        }

        static bool IsRegistered(const std::string &name) {
            return Registry().find(name) != Registry().end();
        }

        static void Display() {
            fmt::print("Registered {}:\n", Base::Name());
            for(const auto &registered : Registry())
                fmt::print("  - {}\n", registered.first);
        }
};

template<typename Base, typename Derived, typename ...Args>
class Registrable {
    protected:
        Registrable() = default;
        virtual ~Registrable() {
            if(!m_registered)
                spdlog::error("Error registering");
        }

        static bool Register() {
            Factory<Base, Args...>::template Register<Derived>(Derived::Name());
            return true;
        }

    private:
        static const bool m_registered;
};
template<typename Base, typename Derived, typename ...Args>
const bool Registrable<Base, Derived, Args...>::m_registered = Registrable<Base, Derived, Args...>::Register();


}

#endif
