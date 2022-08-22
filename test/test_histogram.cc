#include "catch2/catch.hpp"

#include "fmt/format.h"
#include "Achilles/Histogram.hh"

#include <sstream>

TEST_CASE("Histogram Stats", "[Histogram]") {
    achilles::Histogram hist(11, -0.5, 10.5, "test");
    achilles::Histogram hist2({-0.5, 0.5, 1.5, 2.5,
                              3.5, 4.5, 5.5, 6.5,
                              7.5, 8.5, 9.5, 10.5}, "test2");

    hist.Fill(0, 10);
    hist.Fill(1, 9);
    hist.Fill(2, 8);
    hist.Fill(3, 7);
    hist.Fill(4, 6);
    hist.Fill(5, 5);
    hist.Fill(6, 4);
    hist.Fill(7, 3);
    hist.Fill(8, 2);
    hist.Fill(9, 1);
    hist.Fill(10, 0);
    hist2.Fill(0, 10);
    hist2.Fill(1, 9);
    hist2.Fill(2, 8);
    hist2.Fill(3, 7);
    hist2.Fill(4, 6);
    hist2.Fill(5, 5);
    hist2.Fill(6, 4);
    hist2.Fill(7, 3);
    hist2.Fill(8, 2);
    hist2.Fill(9, 1);
    hist2.Fill(10, 0);

    SECTION("Test Integral") {
        CHECK(hist.Integral() == 55);
        CHECK(hist.Integral(0, 5) == 45);

        CHECK(hist.Integral() == -hist.Integral(11, 0));
        CHECK_THROWS_WITH(hist.Integral(12, 2), "Invalid range for histogram integration");
        CHECK_THROWS_WITH(hist.Integral(0, 12), "Invalid range for histogram integration");
    }

    SECTION("Test Normalize and Scale") {
        hist.Scale(1.0/hist.Integral());
        CHECK(hist.Integral() == 1.0);
        hist2.Normalize();
        CHECK(hist2.Integral() == 1.0);
    }
}
