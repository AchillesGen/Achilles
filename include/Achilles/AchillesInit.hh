#ifndef ACHILLES_INIT_HH
#define ACHILLES_INIT_HH

#include <filesystem>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace achilles {

void GenerateEvents(const std::string &, std::vector<std::string> &, int, int, const std::string &,
                    bool);

void InitializePaths(const fs::path &root, const fs::path &libs, const fs::path &data,
                     const fs::path &flux);
} // namespace achilles

#endif // ACHILLES_INIT_HH
