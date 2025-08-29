#ifndef PLUGIN_MANAGER
#define PLUGIN_MANAGER

#include <array>
#include <string>
#include <vector>

namespace achilles {

namespace Plugin {

class Manager {
  public:
    Manager() { Open(); }
    ~Manager() { Close(); }

  private:
    void Open();
    void Close();
    bool ValidateVersion(const std::array<int, 3> &) const;
    std::vector<void *> m_handles;
};

} // namespace Plugin

} // namespace achilles

#endif
