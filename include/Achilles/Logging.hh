// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef LOGGING_HH
#define LOGGING_HH

#include "Achilles/Console.hh"
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
    // Drop the progress bar so the messages, and the shell prompt after them, start on a
    // clean line instead of after the bar and the echoed ^C
    achilles::console::LiveLine::Instance().Clear();
    spdlog::error("Achilles has recieved a {} signal", signal_str);
    spdlog::error("Achilles attempting to exit gracefully!");
    spdlog::get("achilles")->flush();
    exit(-1);
}

inline void CreateLogger(int level, int log_level, int flush_time,
                         std::string destination = "achilles.log") {
    achilles::console::SetupLoggers(level, log_level, flush_time, destination);
}

#endif
