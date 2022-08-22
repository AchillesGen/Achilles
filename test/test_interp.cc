#include "catch2/catch.hpp"

#include "fmt/ostream.h"
#include "Achilles/Interpolation.hh"
#include "Achilles/Utilities.hh"

double TestFunc(const double &x, size_t power) {
    return pow(x, static_cast<double>(power));
}

double TestFunc(const double &x, const double &y, size_t power) {
    return pow(x+y, static_cast<double>(power));
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
        CHECK_THROWS_WITH(achilles::Interp1D(x3, y), "Input and output arrays must be the same size.");
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
        CHECK_THROWS_WITH(achilles::Interp2D(x3, y, z), "Input and output arrays must be the same size.");
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
            auto idxHigh = static_cast<size_t>(std::distance(x.begin(),
                        std::upper_bound(x.begin(), x.end(), xi)));
            auto idxLow = idxHigh-1;
            double val = xi - x[idxLow] < x[idxHigh] - xi ? x[idxLow] : x[idxHigh];
            CHECK(interp(xi) == Approx(TestFunc(val, 1)).epsilon(acc)); 
        }
    }

    SECTION("Linear") {
        std::vector<double> y;
        for(const auto &xi : x) y.emplace_back(TestFunc(xi, 1));
        achilles::Interp1D interp(x, y, achilles::InterpolationType::Polynomial);
        interp.SetPolyOrder(1);

        for(const auto &xi : x_) {
            CHECK(interp(xi) == Approx(TestFunc(xi, 1)).epsilon(acc)); 
        }
    }

    SECTION("Quadratic") {
        std::vector<double> y;
        for(const auto &xi : x) y.emplace_back(TestFunc(xi, 2));
        achilles::Interp1D interp(x, y, achilles::InterpolationType::Polynomial);
        interp.SetPolyOrder(2);

        for(const auto &xi : x_) {
            CHECK(interp(xi) == Approx(TestFunc(xi, 2)).epsilon(acc)); 
        }
    }

    SECTION("Cubic") {
        std::vector<double> y;
        for(const auto &xi : x) y.emplace_back(TestFunc(xi, 3));
        achilles::Interp1D interp(x, y, achilles::InterpolationType::Polynomial);
        interp.SetPolyOrder(3);

        for(const auto &xi : x_) {
            CHECK(interp(xi) == Approx(TestFunc(xi, 3)).epsilon(acc)); 
        }
    }

    SECTION("Cubic Spline") {
        std::vector<double> y;
        for(const auto &xi : x) y.emplace_back(TestFunc(xi, 3));
        achilles::Interp1D interp(x, y, achilles::InterpolationType::CubicSpline);
        CHECK_THROWS_WITH(interp(1), "Interpolation is not initialized!");
        interp.CubicSpline();

        for(const auto &xi : x_) {
            CHECK(interp(xi) == Approx(TestFunc(xi, 3)).epsilon(acc)); 
        }
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
            for(const auto &yi : y)
                z.emplace_back(TestFunc(xi, yi, 1));
        achilles::Interp2D interp(x, y, z, achilles::InterpolationType::NearestNeighbor);

        for(const auto &xi : x_) {
            auto idxHighX = static_cast<size_t>(std::distance(x.begin(),
                        std::upper_bound(x.begin(), x.end(), xi)));
            auto idxLowX = idxHighX-1;
            double valx = xi - x[idxLowX] < x[idxHighX] - xi ? x[idxLowX] : x[idxHighX];
            for(const auto &yi : y_) {
                auto idxHighY = static_cast<size_t>(std::distance(y.begin(),
                            std::upper_bound(y.begin(), y.end(), yi)));
                auto idxLowY = idxHighY-1;
                double valy = yi - y[idxLowY] < y[idxHighY] - yi ? y[idxLowY] : y[idxHighY];
                CHECK(interp(xi, yi) == Approx(TestFunc(valx, valy, 1)).epsilon(acc)); 
            }
        }
    }

    SECTION("Bilinear") {
        std::vector<double> z;
        for(const auto &xi : x) 
            for(const auto &yi : y)
                z.emplace_back(TestFunc(xi, yi, 1));
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
            for(const auto &yi : y)
                z.emplace_back(TestFunc(xi, yi, 2));
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
            for(const auto &yi : y)
                z.emplace_back(TestFunc(xi, yi, 3));
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
            for(const auto &yi : y)
                z.emplace_back(TestFunc(xi, yi, 3));
        achilles::Interp2D interp(x, y, z, achilles::InterpolationType::CubicSpline);
        interp.BicubicSpline();

        for(const auto &xi : x_) {
            for(const auto &yi : y_) {
                CHECK(interp(xi, yi) == Approx(TestFunc(xi, yi, 3)).epsilon(acc)); 
            }
        }

    }
}
