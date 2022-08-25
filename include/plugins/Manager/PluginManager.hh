#ifndef PLUGIN_MANAGER
#define PLUGIN_MANAGER

#include <string>
#include <vector>
#include "Achilles/Utilities.hh"
#include "Achilles/System.hh"

namespace achilles {

namespace Plugin {

class Manager {
    public:
        Manager() { Open(); }
        ~Manager() { Close(); }

    private:
        std::vector<std::string> GetPluginPaths();
        void Open();
        void Close();
        std::vector<void*> m_handles;
};

}

}

#endif
