#ifndef LOGGING_HH
#define LOGGING_HH

#include "spdlog/spdlog.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/sinks/stdout_color_sinks.h"

void CreateLogger(int level, int flush_time) {
    auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
    console_sink->set_level(static_cast<spdlog::level::level_enum>(level));

    auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>("nuchic.log", true);
    file_sink->set_level(spdlog::level::trace);

    spdlog::sinks_init_list sink_list = {file_sink, console_sink};
    auto logger = std::make_shared<spdlog::logger>("nuchic", sink_list);
    logger -> flush_on(spdlog::level::warn);
    spdlog::set_default_logger(logger);
    spdlog::flush_every(std::chrono::seconds(flush_time));
}

#endif
