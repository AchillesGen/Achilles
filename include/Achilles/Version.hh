#ifndef VERSION_HH
#define VERSION_HH

#define ACHILLES_VERSION "1.0.0"
#define ACHILLES_VERSION_MAJOR 1
#define ACHILLES_VERSION_MINOR 0
#define ACHILLES_VERSION_PATCH 0

#include <array>

static constexpr std::array<int, 3> CurrentVersion{ACHILLES_VERSION_MAJOR, ACHILLES_VERSION_MINOR,
                                                   ACHILLES_VERSION_PATCH};

#define EXPORT_ACHILLES_VERSION()                     \
    extern "C" void ExpectedVersion(int version[3]) { \
        version[0] = ACHILLES_VERSION_MAJOR;          \
        version[1] = ACHILLES_VERSION_MINOR;          \
        version[2] = ACHILLES_VERSION_PATCH;          \
    }

#endif
