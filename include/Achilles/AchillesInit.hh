// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef ACHILLES_INIT_HH
#define ACHILLES_INIT_HH

#include <filesystem>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace achilles {

/** Python: Achilles.initialize_logging(verbosity,log_verbosity,logfile) */
void InitializeLogging(int, int, std::string);

/** Python: Achilles.generate_events(runcard,sherpa_args,batch_mode) */
void GenerateEvents(const std::string &, std::vector<std::string> &, bool);

/** Python: Achilles.set_paths(root,libs,data,flux) */
void InitializePaths(const fs::path &root, const fs::path &libs, const fs::path &data,
                     const fs::path &flux);

} // namespace achilles

#endif // ACHILLES_INIT_HH
