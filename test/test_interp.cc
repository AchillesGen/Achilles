#include "catch2/catch.hpp"

#include "Achilles/Interpolation.hh"
#include "Achilles/Utilities.hh"
#include "fmt/ostream.h"
#include <iostream>
#include <sstream>

double TestFunc(const double &x, size_t power) {
    return pow(x, static_cast<double>(power));
}

double TestFuncInt(double min, double max, size_t power) {
    double exponent = static_cast<double>(power + 1);
    return pow(max, exponent) / (exponent)-pow(min, exponent) / (exponent);
}

double TestFunc(const double &x, const double &y, size_t power) {
    return pow(x + y, static_cast<double>(power));
}

constexpr double acc = 1e-4;

TEST_CASE("Errors 1D", "[Interp]") {
    SECTION("Valid Input") {
        const std::vector<double> x = {1, 0, 2, -2};
        const std::vector<double> x2 = {1, 1};
        const std::vector<double> x3 = {0, 1};
        const std::vector<double> y = {0, 1, 2};

        CHECK_THROWS_WITH(achilles::Interp1D(x, y), "Inputs must be increasing.");
        CHECK_THROWS_WITH(achilles::Interp1D(x2, y), "Inputs must all be unique.");
        CHECK_THROWS_WITH(achilles::Interp1D(x3, y),
                          "Input and output arrays must be the same size.");
    }

    SECTION("No extrapolation allowed") {
        const std::vector<double> x = {0, 1, 2};
        const std::vector<double> y = {0, 1, 2};

        achilles::Interp1D interp(x, y, achilles::InterpolationType::Polynomial);
        interp.SetPolyOrder(2);

        CHECK_THROWS_WITH(interp(-1),
                          fmt::format("Input ({}) less than minimum value ({})", -1, 0));
        CHECK_THROWS_WITH(interp(3),
                          fmt::format("Input ({}) greater than maximum value ({})", 3, 2));
    }
}

TEST_CASE("Errors 2D", "[Interp]") {
    SECTION("Valid Input") {
        const std::vector<double> x = {1, 0};
        const std::vector<double> x2 = {0, 0};
        const std::vector<double> x3 = {0, 1};
        const std::vector<double> y = {0, 1, 2};
        const std::vector<double> y2 = {0, 1, 0};
        const std::vector<double> y3 = {0, 0};
        const std::vector<double> z = {0, 1, 2};

        CHECK_THROWS_WITH(achilles::Interp2D(x, y, z), "Inputs must be increasing.");
        CHECK_THROWS_WITH(achilles::Interp2D(x2, y, z), "Inputs must all be unique.");
        CHECK_THROWS_WITH(achilles::Interp2D(x3, y2, z), "Inputs must be increasing.");
        CHECK_THROWS_WITH(achilles::Interp2D(x3, y3, z), "Inputs must all be unique.");
        CHECK_THROWS_WITH(achilles::Interp2D(x3, y, z),
                          "Input and output arrays must be the same size.");
    }

    SECTION("No extrapolation allowed") {
        const std::vector<double> x = {0, 1, 2};
        const std::vector<double> y = {0, 1, 2};
        const std::vector<double> z = {0, 1, 2, 3, 4, 5, 6, 7, 8};

        achilles::Interp2D interp(x, y, z, achilles::InterpolationType::Polynomial);
        interp.SetPolyOrder(2, 2);

        CHECK_THROWS_WITH(interp(-1, 0),
                          fmt::format("Input ({}) less than minimum x value ({})", -1, 0));
        CHECK_THROWS_WITH(interp(3, 0),
                          fmt::format("Input ({}) greater than maximum x value ({})", 3, 2));
        CHECK_THROWS_WITH(interp(0, -1),
                          fmt::format("Input ({}) less than minimum y value ({})", -1, 0));
        CHECK_THROWS_WITH(interp(0, 3),
                          fmt::format("Input ({}) greater than maximum y value ({})", 3, 2));
    }
}

TEST_CASE("One Dimensional", "[Interp]") {
    const std::vector<double> x = achilles::Linspace(0, 10, 97);
    const std::vector<double> x_ = achilles::Linspace(0, 10, 101);

    SECTION("Nearest Neighbor") {
        std::vector<double> y;
        for(const auto &xi : x) y.emplace_back(TestFunc(xi, 1));
        achilles::Interp1D interp(x, y, achilles::InterpolationType::NearestNeighbor);

        for(const auto &xi : x_) {
            auto idxHigh = static_cast<size_t>(
                std::distance(x.begin(), std::upper_bound(x.begin(), x.end(), xi)));
            auto idxLow = idxHigh - 1;
            double val = xi - x[idxLow] < x[idxHigh] - xi ? x[idxLow] : x[idxHigh];
            CHECK(interp(xi) == Approx(TestFunc(val, 1)).epsilon(acc));
        }

        CHECK_THAT(interp.Integrate(), Catch::Matchers::WithinRel(TestFuncInt(0, 10, 1), 1e-4));
    }

    SECTION("Linear") {
        std::vector<double> y;
        for(const auto &xi : x) y.emplace_back(TestFunc(xi, 1));
        achilles::Interp1D interp(x, y, achilles::InterpolationType::Polynomial);
        interp.SetPolyOrder(1);

        for(const auto &xi : x_) { CHECK(interp(xi) == Approx(TestFunc(xi, 1)).epsilon(acc)); }
        CHECK_THAT(interp.Integrate(), Catch::Matchers::WithinRel(TestFuncInt(0, 10, 1), 1e-4));
    }

    SECTION("Quadratic") {
        std::vector<double> y;
        for(const auto &xi : x) y.emplace_back(TestFunc(xi, 2));
        achilles::Interp1D interp(x, y, achilles::InterpolationType::Polynomial);
        interp.SetPolyOrder(2);

        for(const auto &xi : x_) { CHECK(interp(xi) == Approx(TestFunc(xi, 2)).epsilon(acc)); }
        CHECK_THAT(interp.Integrate(), Catch::Matchers::WithinRel(TestFuncInt(0, 10, 2), 1e-4));
    }

    SECTION("Cubic") {
        std::vector<double> y;
        for(const auto &xi : x) y.emplace_back(TestFunc(xi, 3));
        achilles::Interp1D interp(x, y, achilles::InterpolationType::Polynomial);
        interp.SetPolyOrder(3);

        for(const auto &xi : x_) { CHECK(interp(xi) == Approx(TestFunc(xi, 3)).epsilon(acc)); }
    }

    SECTION("Cubic Spline") {
        std::vector<double> y;
        for(const auto &xi : x) y.emplace_back(TestFunc(xi, 3));
        achilles::Interp1D interp(x, y, achilles::InterpolationType::CubicSpline);
        CHECK_THROWS_WITH(interp(1), "Interpolation is not initialized!");
        interp.CubicSpline();

        for(const auto &xi : x_) { CHECK(interp(xi) == Approx(TestFunc(xi, 3)).epsilon(acc)); }
        CHECK_THAT(interp.Integrate(), Catch::Matchers::WithinRel(TestFuncInt(0, 10, 3), 1e-4));
    }
}

TEST_CASE("Derivative and extrapolation", "[Interp]") {
    // A straight line y = 2x + 1: the cubic spline is exact, its derivative is 2
    // everywhere, and linear extrapolation continues the line.
    const std::vector<double> x = achilles::Linspace(0, 10, 51);
    std::vector<double> y;
    for(const auto &xi : x) y.emplace_back(2 * xi + 1);
    achilles::Interp1D interp(x, y, achilles::InterpolationType::CubicSpline);
    interp.CubicSpline();

    SECTION("derivative on a line is the slope") {
        CHECK(interp.Derivative(3.0) == Approx(2.0).margin(1e-6));
        CHECK(interp.Derivative(0.0) == Approx(2.0).margin(1e-6));  // left end
        CHECK(interp.Derivative(10.0) == Approx(2.0).margin(1e-6)); // right end
    }

    SECTION("extrapolation is opt-in") {
        CHECK_THROWS_AS(interp(-1.0), std::domain_error);
        CHECK_THROWS_AS(interp(11.0), std::domain_error);
    }

    SECTION("linear extrapolation follows the boundary slope") {
        interp.AllowExtrapolation(true);
        CHECK(interp(-2.0) == Approx(2 * -2.0 + 1).margin(1e-6));
        CHECK(interp(13.0) == Approx(2 * 13.0 + 1).margin(1e-6));
    }
}

TEST_CASE("Monotone cubic", "[Interp]") {
    SECTION("exact on a straight line") {
        const std::vector<double> x = achilles::Linspace(0, 10, 11);
        std::vector<double> y;
        for(const auto &xi : x) y.emplace_back(2 * xi + 1);
        achilles::Interp1D interp(x, y);
        interp.MonotoneSpline();
        for(double xq = 0; xq <= 10; xq += 0.13)
            CHECK(interp(xq) == Approx(2 * xq + 1).margin(1e-9));
        CHECK(interp.Derivative(4.2) == Approx(2.0).margin(1e-9));
    }

    SECTION("never overshoots its own knots") {
        // A step is what a natural cubic rings against: it undershoots below the
        // low plateau and overshoots above the high one. The monotone form
        // cannot leave the [min, max] of the knots it interpolates between.
        const std::vector<double> x = achilles::Linspace(0, 10, 11);
        const std::vector<double> y = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1};
        achilles::Interp1D mono(x, y);
        mono.MonotoneSpline();
        achilles::Interp1D natural(x, y);
        natural.CubicSpline();

        bool natural_overshoots = false;
        for(double xq = 0; xq <= 10; xq += 0.01) {
            CHECK(mono(xq) >= 0.0);
            CHECK(mono(xq) <= 1.0);
            if(natural(xq) < -1e-6 || natural(xq) > 1 + 1e-6) natural_overshoots = true;
        }
        CHECK(natural_overshoots);
    }

    SECTION("preserves monotonicity of noisy but rising data") {
        const std::vector<double> x = {0, 1, 2, 3, 4, 5, 6};
        const std::vector<double> y = {0, 1, 1.05, 3, 3.02, 5, 9};
        achilles::Interp1D interp(x, y);
        interp.MonotoneSpline();
        double prev = interp(0.0);
        for(double xq = 0.005; xq <= 6; xq += 0.005) {
            CHECK(interp(xq) >= prev - 1e-12);
            prev = interp(xq);
        }
    }

    SECTION("survives a cache round trip") {
        const std::vector<double> x = {0, 1, 2, 3, 4};
        const std::vector<double> y = {0, 2, 2.1, 5, 9};
        achilles::Interp1D interp(x, y);
        interp.MonotoneSpline();

        std::stringstream ss;
        interp.SaveState(ss);
        achilles::Interp1D loaded;
        loaded.LoadState(ss);
        for(double xq = 0; xq <= 4; xq += 0.05)
            CHECK(loaded(xq) == Approx(interp(xq)).margin(1e-12));
    }
}

TEST_CASE("Two Dimensional", "[Interp]") {
    const std::vector<double> x = achilles::Linspace(0, 11, 97);
    const std::vector<double> y = achilles::Linspace(0, 11, 97);
    const std::vector<double> x_ = achilles::Linspace(0, 10, 11);
    const std::vector<double> y_ = achilles::Linspace(0, 10, 11);

    SECTION("Nearest Neighbor") {
        std::vector<double> z;
        for(const auto &xi : x)
            for(const auto &yi : y) z.emplace_back(TestFunc(xi, yi, 1));
        achilles::Interp2D interp(x, y, z, achilles::InterpolationType::NearestNeighbor);

        for(const auto &xi : x_) {
            auto idxHighX = static_cast<size_t>(
                std::distance(x.begin(), std::upper_bound(x.begin(), x.end(), xi)));
            auto idxLowX = idxHighX - 1;
            double valx = xi - x[idxLowX] < x[idxHighX] - xi ? x[idxLowX] : x[idxHighX];
            for(const auto &yi : y_) {
                auto idxHighY = static_cast<size_t>(
                    std::distance(y.begin(), std::upper_bound(y.begin(), y.end(), yi)));
                auto idxLowY = idxHighY - 1;
                double valy = yi - y[idxLowY] < y[idxHighY] - yi ? y[idxLowY] : y[idxHighY];
                CHECK(interp(xi, yi) == Approx(TestFunc(valx, valy, 1)).epsilon(acc));
            }
        }
    }

    SECTION("Bilinear") {
        std::vector<double> z;
        for(const auto &xi : x)
            for(const auto &yi : y) z.emplace_back(TestFunc(xi, yi, 1));
        achilles::Interp2D interp(x, y, z, achilles::InterpolationType::Polynomial);
        interp.SetPolyOrder(1, 1);

        for(const auto &xi : x_) {
            for(const auto &yi : y_) {
                CHECK(interp(xi, yi) == Approx(TestFunc(xi, yi, 1)).epsilon(acc));
            }
        }
    }

    SECTION("Biquadratic") {
        std::vector<double> z;
        for(const auto &xi : x)
            for(const auto &yi : y) z.emplace_back(TestFunc(xi, yi, 2));
        achilles::Interp2D interp(x, y, z, achilles::InterpolationType::Polynomial);
        interp.SetPolyOrder(2, 2);

        for(const auto &xi : x_) {
            for(const auto &yi : y_) {
                CHECK(interp(xi, yi) == Approx(TestFunc(xi, yi, 2)).epsilon(acc));
            }
        }
    }

    SECTION("Bicubic") {
        std::vector<double> z;
        for(const auto &xi : x)
            for(const auto &yi : y) z.emplace_back(TestFunc(xi, yi, 3));
        achilles::Interp2D interp(x, y, z, achilles::InterpolationType::Polynomial);
        interp.SetPolyOrder(3, 3);

        for(const auto &xi : x_) {
            for(const auto &yi : y_) {
                CHECK(interp(xi, yi) == Approx(TestFunc(xi, yi, 3)).epsilon(acc));
            }
        }
    }

    SECTION("Bicubic Spline") {
        std::vector<double> z;
        for(const auto &xi : x)
            for(const auto &yi : y) z.emplace_back(TestFunc(xi, yi, 3));
        achilles::Interp2D interp(x, y, z, achilles::InterpolationType::CubicSpline);
        interp.BicubicSpline();

        for(const auto &xi : x_) {
            for(const auto &yi : y_) {
                CHECK(interp(xi, yi) == Approx(TestFunc(xi, yi, 3)).epsilon(acc));
            }
        }
    }
}
