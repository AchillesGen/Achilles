#include "spdlog/spdlog.h"

extern "C" {

    void LogDebug(const char *msg) {
        spdlog::debug(msg); 
    }

    void LogInfo(const char *msg) {
        spdlog::info(msg);
    }

    void LogError(const char *msg) {
        spdlog::error(msg);
    }

    void LogWarn(const char *msg) {
        spdlog::warn(msg);
    }

    void LogCritical(const char *msg) {
        spdlog::critical(msg);
    }

}
