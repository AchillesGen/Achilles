#include "pybind11/pybind11.h"

#include "spdlog/spdlog.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/sinks/stdout_color_sinks.h"

namespace py = pybind11;

PYBIND11_MODULE(logger, m) {
    m.def("init",  [](const std::string &logname, const std::string &logfile) {
        auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
        console_sink->set_level(spdlog::level::info);
        console_sink->set_pattern("[%^%l%$] %v");

        auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(logfile, true); 
        file_sink->set_level(spdlog::level::debug);
    
        spdlog::logger logger(logname, {console_sink, file_sink});
        auto logger_ptr = std::make_shared<spdlog::logger>(logger);
        logger_ptr -> set_level(spdlog::level::debug);
        spdlog::set_default_logger(logger_ptr);
    });

    py::enum_<spdlog::level::level_enum>(m, "level")
        .value("trace", spdlog::level::level_enum::trace)
        .value("debug", spdlog::level::level_enum::debug)
        .value("info", spdlog::level::level_enum::info)
        .value("warning", spdlog::level::level_enum::warn)
        .value("error", spdlog::level::level_enum::err)
        .value("critical", spdlog::level::level_enum::critical)
        .value("off", spdlog::level::level_enum::off);

    m.def("set_level", [](const spdlog::level::level_enum& level) {
        spdlog::set_level(level);
    });
    m.def("set_pattern", [](const std::string& pattern) { spdlog::set_pattern(pattern); });

    m.def("log", [](const spdlog::level::level_enum& level, const std::string& message) {
        spdlog::log(level, message);
    });

    m.def("trace", [](const std::string &message) { spdlog::trace(message); });
    m.def("debug", [](const std::string &message) { spdlog::debug(message); });
    m.def("info", [](const std::string &message) { spdlog::info(message); });
    m.def("warning", [](const std::string &message) { spdlog::warn(message); });
    m.def("error", [](const std::string &message) { spdlog::error(message); });
    m.def("critical", [](const std::string &message) { spdlog::critical(message); });
}
