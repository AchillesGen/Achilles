#include "catch2/catch_approx.hpp"
#include "catch2/catch_test_macros.hpp"
#include "catch2/generators/catch_generators_adapters.hpp"
#include "catch2/generators/catch_generators_random.hpp"

#include "Achilles/Autodiff.hh"
#include "catch_utils.hh"

TEST_CASE("Dual Numbers", "[Utilities]") {
    achilles::Dual x(GENERATE(take(10, random<double>(0, 100))));

    SECTION("Operators") {
        achilles::Dual y(GENERATE(take(10, random<double>(-100, 100))));
        double a(GENERATE(take(10, random<double>(-100, 100))));

        auto z = +x;
        CHECK(z.Value() == x.Value());
        CHECK(z.Derivative() == x.Derivative());

        z = -x;
        CHECK(z.Value() == -x.Value());
        CHECK(z.Derivative() == -x.Derivative());

        z = x + y;
        CHECK(z.Value() == x.Value() + y.Value());
        CHECK(z.Derivative() == Catch::Approx(2));

        z = x - y;
        CHECK(z.Value() == x.Value() - y.Value());
        CHECK(z.Derivative() == Catch::Approx(0));

        z = x * y;
        CHECK(z.Value() == x.Value() * y.Value());
        CHECK(z.Derivative() == Catch::Approx(y.Value() + x.Value()));

        z = x / y;
        CHECK(z.Value() == x.Value() / y.Value());
        CHECK(z.Derivative() == Catch::Approx((y.Value() - x.Value()) / y.Value() / y.Value()));

        z = a + x;
        CHECK(z.Value() == x.Value() + a);
        CHECK(z.Derivative() == Catch::Approx(1));

        z = x + a;
        CHECK(z.Value() == x.Value() + a);
        CHECK(z.Derivative() == Catch::Approx(1));

        z = x - a;
        CHECK(z.Value() == x.Value() - a);
        CHECK(z.Derivative() == Catch::Approx(1));

        z = a - x;
        CHECK(z.Value() == a - x.Value());
        CHECK(z.Derivative() == Catch::Approx(-1));

        z = x * a;
        auto z2 = a * x;
        CHECK(z.Value() == z2.Value());
        CHECK(z.Derivative() == Catch::Approx(z2.Derivative()));
        CHECK(z.Derivative() == Catch::Approx(a));

        z = x / a;
        CHECK(z.Value() == x.Value() / a);
        CHECK(z.Derivative() == Catch::Approx(1.0 / a));

        z = a / x;
        CHECK(z.Value() == a / x.Value());
        CHECK(z.Derivative() == Catch::Approx(-a / x.Value() / x.Value()));
    }

    SECTION("Special functions") {
        auto z = sin(x);
        CHECK(z.Value() == sin(x.Value()));
        CHECK(z.Derivative() == Catch::Approx(cos(x.Value())));

        z = cos(x);
        CHECK(z.Value() == cos(x.Value()));
        CHECK(z.Derivative() == Catch::Approx(-sin(x.Value())));

        z = tan(x);
        CHECK(z.Value() == tan(x.Value()));
        CHECK(z.Derivative() == Catch::Approx(1.0 / std::pow(std::cos(x.Value()), 2)));

        z = exp(x);
        CHECK(z.Value() == exp(x.Value()));
        CHECK(z.Derivative() == Catch::Approx(std::exp(x.Value())));

        z = log(x);
        CHECK(z.Value() == log(x.Value()));
        CHECK(z.Derivative() == Catch::Approx(1 / x.Value()));

        z = abs(x);
        CHECK(z.Value() == std::abs(x.Value()));
        const double sign = x.Value() == 0 ? 0 : x.Value() / std::abs(x.Value());
        CHECK(z.Derivative() == Catch::Approx(sign));

        z = cosh(x);
        CHECK(z.Value() == cosh(x.Value()));
        CHECK(z.Derivative() == Catch::Approx(std::sinh(x.Value())));

        z = sinh(x);
        CHECK(z.Value() == sinh(x.Value()));
        CHECK(z.Derivative() == Catch::Approx(std::cosh(x.Value())));

        z = sech(x);
        CHECK(z.Value() == 1.0 / cosh(x.Value()));
        CHECK(z.Derivative() == Catch::Approx(-std::tanh(x.Value()) / std::cosh(x.Value())));

        z = pow(x, 3);
        CHECK(z.Value() == pow(x.Value(), 3));
        CHECK(z.Derivative() == Catch::Approx(3 * pow(x.Value(), 2)));
    }

    SECTION("Chain rule") {
        double xval = x.Value();
        auto z = sin(abs(log(x))) * pow(x, 3) / exp(x) + x;
        CHECK(z.Value() ==
              sin(std::abs(log(x.Value()))) * pow(x.Value(), 3) / exp(x.Value()) + x.Value());
        CHECK(z.Derivative() ==
              Catch::Approx(-exp(-xval) * pow(xval, 3) * sin(std::abs(log(xval))) +
                            3 * exp(-xval) * pow(xval, 2) * sin(std::abs(log(xval))) +
                            exp(-xval) * pow(xval, 2) * log(xval) * cos(std::abs(log(xval))) /
                                std::abs(log(xval)) +
                            1));
    }
}
