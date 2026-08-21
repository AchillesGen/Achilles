#include "catch2/catch.hpp"

#include "Achilles/Profiler.hh"

#include <ctime>
#include <thread>

TEST_CASE("Profiler accumulates sections", "[Profiler]") {
    auto &profiler = achilles::Profiler::Instance();
    profiler.Reset();

    SECTION("Repeated sections are accumulated") {
        profiler.Record("section", 1.5);
        profiler.Record("section", 2.5);
        profiler.Record("other", 1.0);

        CHECK(profiler.Time("section") == Approx(4.0));
        CHECK(profiler.Calls("section") == 2);
        CHECK(profiler.Time("other") == Approx(1.0));
        CHECK(profiler.Calls("other") == 1);
    }

    SECTION("Sections keep their insertion order") {
        profiler.Record("second", 1.0);
        profiler.Record("first", 1.0);
        profiler.Record("second", 1.0);

        CHECK(profiler.Sections().at("second").order == 0);
        CHECK(profiler.Sections().at("first").order == 1);
    }

    SECTION("Unknown sections have no time") {
        CHECK(profiler.Time("missing") == Approx(0.0));
        CHECK(profiler.Calls("missing") == 0);
    }

    SECTION("ScopedTimer records the enclosing scope") {
        const auto sleep = std::chrono::milliseconds(2);
        {
            achilles::ScopedTimer timer{"scope"};
            std::this_thread::sleep_for(sleep);
        }
        CHECK(profiler.Calls("scope") == 1);
        CHECK(profiler.Time("scope") >= std::chrono::duration<double>(sleep).count());
    }

    SECTION("ScopedTimer only records once") {
        {
            achilles::ScopedTimer timer{"scope"};
            timer.Stop();
            timer.Stop();
        }
        CHECK(profiler.Calls("scope") == 1);
    }

    profiler.Reset();
}

TEST_CASE("Profiler formats timestamps", "[Profiler]") {
    // 2026-08-21 11:09:47 local time
    std::tm tm{};
    tm.tm_year = 126;
    tm.tm_mon = 7;
    tm.tm_mday = 21;
    tm.tm_hour = 11;
    tm.tm_min = 9;
    tm.tm_sec = 47;
    tm.tm_isdst = -1;
    const auto time = std::chrono::system_clock::from_time_t(std::mktime(&tm));

    CHECK(achilles::Profiler::FormatTimestamp(time) == "Fri Aug 21 11:09:47 2026");
}

TEST_CASE("Profiler formats times", "[Profiler]") {
    CHECK(achilles::Profiler::FormatTime(1.5e-5) == "15.0 us");
    CHECK(achilles::Profiler::FormatTime(0.25) == "250.0 ms");
    CHECK(achilles::Profiler::FormatTime(12.5) == "12.500 s");
    CHECK(achilles::Profiler::FormatTime(75) == "1m 15s");
    CHECK(achilles::Profiler::FormatTime(3725) == "1h 02m 05s");
}
