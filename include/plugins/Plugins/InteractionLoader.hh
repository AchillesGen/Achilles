#ifndef INTERACTIONLOADER_HH
#define INTERACTIONLOADER_HH

#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

#if __cplusplus == 201402L
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#elif __cplusplus >= 201703L
#include <filesystem>
namespace fs = std::filesystem;
#else
#error Plugin Library: Currently only works with c++14 or newer
#endif

namespace achilles {

class Interactions;
class InteractionFactory;

class InteractionLoader {
    public:
        InteractionLoader() = delete;
        static bool LoadInteractions(const std::vector<std::string>&);

    private:
        void _registerInteraction();
        void _loadInteractionPlugins();

        static bool m_init;
        static std::vector<fs::path> plugins;
        static std::set<Interactions*> loadedPlugins;
};

} // end achilles namespace

#endif
