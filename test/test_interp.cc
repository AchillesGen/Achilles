#include "catch2/catch.hpp"

#include "fmt/ostream.h"
#include "nuchic/Interpolation.hh"
#include "nuchic/Utilities.hh"

double TestFunc(const double &x, size_t power) {
    return pow(x, static_cast<double>(power));
}

double TestFunc(const double &x, const double &y, size_t power) {
    return pow(x+y, static_cast<double>(power));
}

constexpr double acc = 1e-4;

TEST_CASE("One Dimensional", "[Interp]") {
    const std::vector<double> x = nuchic::Linspace(0, 10, 97);
    const std::vector<double> x_ = nuchic::Linspace(0, 10, 101);

    SECTION("Nearest Neighbor") {
        std::vector<double> y;
        for(const auto &xi : x) y.emplace_back(TestFunc(xi, 1));
        nuchic::Interp1D interp(x, y, nuchic::InterpolationType::NearestNeighbor);

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
        nuchic::Interp1D interp(x, y, nuchic::InterpolationType::Polynomial);
        interp.SetPolyOrder(1);

        for(const auto &xi : x_) {
            CHECK(interp(xi) == Approx(TestFunc(xi, 1)).epsilon(acc)); 
        }
    }

    SECTION("Quadratic") {
        std::vector<double> y;
        for(const auto &xi : x) y.emplace_back(TestFunc(xi, 2));
        nuchic::Interp1D interp(x, y, nuchic::InterpolationType::Polynomial);
        interp.SetPolyOrder(2);

        for(const auto &xi : x_) {
            CHECK(interp(xi) == Approx(TestFunc(xi, 2)).epsilon(acc)); 
        }
    }

    SECTION("Cubic") {
        std::vector<double> y;
        for(const auto &xi : x) y.emplace_back(TestFunc(xi, 3));
        nuchic::Interp1D interp(x, y, nuchic::InterpolationType::Polynomial);
        interp.SetPolyOrder(3);

        for(const auto &xi : x_) {
            CHECK(interp(xi) == Approx(TestFunc(xi, 3)).epsilon(acc)); 
        }
    }

    SECTION("Cubic Spline") {
        std::vector<double> y;
        for(const auto &xi : x) y.emplace_back(TestFunc(xi, 3));
        nuchic::Interp1D interp(x, y, nuchic::InterpolationType::CubicSpline);
        interp.CubicSpline();

        for(const auto &xi : x_) {
            CHECK(interp(xi) == Approx(TestFunc(xi, 3)).epsilon(acc)); 
        }
    }
}

TEST_CASE("Two Dimensional", "[Interp]") {
    const std::vector<double> x = nuchic::Linspace(0, 11, 97);
    const std::vector<double> y = nuchic::Linspace(0, 11, 97);
    const std::vector<double> x_ = nuchic::Linspace(0, 10, 11);
    const std::vector<double> y_ = nuchic::Linspace(0, 10, 11);

    SECTION("Nearest Neighbor") {
        std::vector<double> z;
        for(const auto &xi : x) 
            for(const auto &yi : y)
                z.emplace_back(TestFunc(xi, yi, 1));
        nuchic::Interp2D interp(x, y, z, nuchic::InterpolationType::NearestNeighbor);

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
        nuchic::Interp2D interp(x, y, z, nuchic::InterpolationType::Polynomial);
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
        nuchic::Interp2D interp(x, y, z, nuchic::InterpolationType::Polynomial);
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
        nuchic::Interp2D interp(x, y, z, nuchic::InterpolationType::Polynomial);
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
        nuchic::Interp2D interp(x, y, z, nuchic::InterpolationType::CubicSpline);
        interp.BicubicSpline();

        for(const auto &xi : x_) {
            for(const auto &yi : y_) {
                CHECK(interp(xi, yi) == Approx(TestFunc(xi, yi, 3)).epsilon(acc)); 
            }
        }

    }
}
