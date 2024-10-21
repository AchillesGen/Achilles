#ifndef LOGGING_HH
#define LOGGING_HH

#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/spdlog.h"

#include <csignal>

inline void SignalHandler(int signal) {
    std::string signal_str;
    switch(signal) {
    case SIGTERM:
        signal_str = "SIGTERM";
        break;
    case SIGSEGV:
        signal_str = "SIGSEGV";
        break;
    case SIGINT:
        signal_str = "SIGINT";
        break;
    case SIGABRT:
        signal_str = "SIGABRT";
        break;
    default:
        signal_str = "Unknown";
    }
    spdlog::error("Achilles has recieved a {} signal", signal_str);
    spdlog::error("Achilles attempting to exit gracefully!");
    spdlog::get("achilles")->flush();
    exit(-1);
}

inline void CreateLogger(int level, int flush_time) {
    auto slevel = static_cast<spdlog::level::level_enum>(level);
    auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
    console_sink->set_level(slevel);

    auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>("achilles.log", true);
    file_sink->set_level(slevel);

    spdlog::sinks_init_list sink_list = {file_sink, console_sink};
    auto logger = std::make_shared<spdlog::logger>("achilles", sink_list);
    logger->set_level(slevel);
    logger->flush_on(spdlog::level::warn);
    spdlog::register_logger(logger);
    spdlog::set_default_logger(logger);
    spdlog::flush_every(std::chrono::seconds(flush_time));
    spdlog::set_pattern("[%^%l%$] %v");
}

#endif
