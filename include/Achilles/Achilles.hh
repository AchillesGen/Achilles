#include <string>
#include <vector>

#ifndef ACHILLES_HH
#define ACHILLES_HH

#if TESTING
#define MOCK virtual
#else
#define MOCK
#endif

namespace achilles {
inline void GenerateEvents(const std::string &, std::vector<std::string> &, int, int,
                           const std::string &, bool);
}

#endif
