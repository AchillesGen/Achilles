#include "catch2/catch.hpp"
#include "mock_classes.hh"
#include "Achilles/Potential.hh"
#include "Achilles/Particle.hh"
#include <iostream>

double stencil5(std::function<double(double)> f, double x, double h) {
    return (-f(x + 2*h) + 8*f(x+h) - 8*f(x-h) - f(x-2*h))/(12*h);
}

template<typename T>
T test_derivative(const T& x, const T& y, const T& z) {
    return (x + y + z) * exp(x * y * z);
}

TEST_CASE("CooperPotential::EDAD1 Values", "[Potential]") {
    constexpr double tplab = 100;
    constexpr size_t AA = 12;

    SECTION("Values and Derivatives are consistent") {
        // Results from the fortran code for tplab = 100, r = 0.15, A = 12
        constexpr double rvector = 301.65297587389080, ivector = -84.374377568119442;
        constexpr double rscalar = -384.29852972736001, iscalar = 89.69329129030613;

        auto nucleus = std::make_shared<MockNucleus>();
        REQUIRE_CALL(*nucleus, NNucleons())
            .LR_RETURN((AA))
            .TIMES(AT_LEAST(1));

        achilles::CooperPotential potential;

#ifdef AUTODIFF
        autodiff::real r = 0.15;
        autodiff::real plab = sqrt(pow(tplab + achilles::Constant::mN, 2) - pow(achilles::Constant::mN, 2));

        auto func_rv = [&](const autodiff::real &x, const autodiff::real &y){ return potential.evaluate<autodiff::real>(x, y).rvector; };
        auto func_iv = [&](const autodiff::real &x, const autodiff::real &y){ return potential.evaluate<autodiff::real>(x, y).ivector; };
        auto func_rs = [&](const autodiff::real &x, const autodiff::real &y){ return potential.evaluate<autodiff::real>(x, y).rscalar; };
        auto func_is = [&](const autodiff::real &x, const autodiff::real &y){ return potential.evaluate<autodiff::real>(x, y).iscalar; };
        double drvdp = autodiff::derivative(func_rv, autodiff::wrt(plab), autodiff::at(plab, r));
        double drvdr = autodiff::derivative(func_rv, autodiff::wrt(r), autodiff::at(plab, r));
        double divdp = autodiff::derivative(func_iv, autodiff::wrt(plab), autodiff::at(plab, r));
        double divdr = autodiff::derivative(func_iv, autodiff::wrt(r), autodiff::at(plab, r));
        double drsdp = autodiff::derivative(func_rs, autodiff::wrt(plab), autodiff::at(plab, r));
        double drsdr = autodiff::derivative(func_rs, autodiff::wrt(r), autodiff::at(plab, r));
        double disdp = autodiff::derivative(func_is, autodiff::wrt(plab), autodiff::at(plab, r));
        double disdr = autodiff::derivative(func_is, autodiff::wrt(r), autodiff::at(plab, r));

        std::cout << drvdp << " " << drvdr << "\n";
        std::cout << divdp << " " << divdr << "\n";
        std::cout << drsdp << " " << drsdr << "\n";
        std::cout << disdp << " " << disdr << "\n";

        auto stencilp = potential.derivative_p(nucleus.get(), plab, r, 0.01);
        auto stencilr = potential.derivative_r(nucleus.get(), plab, r, 0.01);

        CHECK(drvdp == Approx(stencilp.rvector));
        CHECK(drvdr == Approx(stencilr.rvector));
        CHECK(divdp == Approx(stencilp.ivector));
        CHECK(divdr == Approx(stencilr.ivector));
        CHECK(drsdp == Approx(stencilp.rscalar));
        CHECK(drsdr == Approx(stencilr.rscalar));
        CHECK(disdp == Approx(stencilp.iscalar));
        CHECK(disdr == Approx(stencilr.iscalar));
#else
        double r = 0.15;
        double plab = sqrt(pow(tplab + achilles::Constant::mN, 2) - pow(achilles::Constant::mN, 2));

#endif
        auto vals = potential(nucleus.get(), plab, r);
        auto stencilp = potential.derivative_p(nucleus.get(), plab, r);
        auto stencilr = potential.derivative_r(nucleus.get(), plab, r);

        // Require to match to 0.01% of Cooper code
        CHECK(vals.rvector == Approx(rvector).epsilon(0.0001));
        CHECK(vals.ivector == Approx(ivector).epsilon(0.0001));
        CHECK(vals.rscalar == Approx(rscalar).epsilon(0.0001));
        CHECK(vals.iscalar == Approx(iscalar).epsilon(0.0001));

        CHECK(stencilp.rvector != 0);
        CHECK(stencilr.rvector != 0);

        std::cout << stencilp.rvector << " " << stencilp.ivector << "\n";
        std::cout << stencilp.rscalar << " " << stencilp.iscalar << "\n";
        std::cout << stencilr.rvector << " " << stencilr.ivector << "\n";
        std::cout << stencilr.rscalar << " " << stencilr.iscalar << "\n";
    }

#if defined(CATCH_CONFIG_ENABLE_BENCHMARKING)

#ifdef AUTODIFF
    BENCHMARK_ADVANCED("Autodiff_test")(Catch::Benchmark::Chronometer meter) {
        autodiff::dual x = 1.0;
        autodiff::dual y = 2.0;
        autodiff::dual z = 3.0;

        autodiff::dual u = test_derivative(x, y, z);
        meter.measure([&]() {
                double dudx = autodiff::derivative(test_derivative<autodiff::dual>, autodiff::wrt(x), autodiff::at(x, y, z));
                double dudy = autodiff::derivative(test_derivative<autodiff::dual>, autodiff::wrt(y), autodiff::at(x, y, z));
                double dudz = autodiff::derivative(test_derivative<autodiff::dual>, autodiff::wrt(z), autodiff::at(x, y, z));
                return u + dudx + dudy + dudz;
        });
    };
#endif

    BENCHMARK_ADVANCED("Stencil_test")(Catch::Benchmark::Chronometer meter) {
        double x = 1.0;
        double y = 2.0;
        double z = 3.0;

        auto func_x = [&](double val) { return test_derivative(val, y, z); };
        auto func_y = [&](double val) { return test_derivative(x, val, z); };
        auto func_z = [&](double val) { return test_derivative(x, y, val); };

        meter.measure([&]() {
                double dudx = stencil5(func_x, x, 0.01);
                double dudy = stencil5(func_y, y, 0.01);
                double dudz = stencil5(func_z, z, 0.01);
                return dudx + dudy + dudz;
        });
    };

#ifdef AUTODIFF
    BENCHMARK_ADVANCED("Autodiff")(Catch::Benchmark::Chronometer meter) {
        autodiff::dual r = 0.15;
        auto nucleus = std::make_shared<MockNucleus>();
        REQUIRE_CALL(*nucleus, NNucleons())
            .LR_RETURN((AA))
            .TIMES(AT_LEAST(1));
        achilles::CooperPotential potential(nucleus);
        autodiff::dual plab = sqrt(pow(tplab + achilles::Constant::mN, 2) - pow(achilles::Constant::mN, 2));

        auto func_rv = [&](const autodiff::dual &x, const autodiff::dual &y){ return potential.evaluate<autodiff::dual>(x, y).rvector; };
        auto func_iv = [&](const autodiff::dual &x, const autodiff::dual &y){ return potential.evaluate<autodiff::dual>(x, y).ivector; };
        auto func_rs = [&](const autodiff::dual &x, const autodiff::dual &y){ return potential.evaluate<autodiff::dual>(x, y).rscalar; };
        auto func_is = [&](const autodiff::dual &x, const autodiff::dual &y){ return potential.evaluate<autodiff::dual>(x, y).iscalar; };

        meter.measure([&]() {
                // potential.evaluate<autodiff::dual>(plab, r);
                autodiff::derivative(func_rv, autodiff::wrt(plab), autodiff::at(plab, r));
                autodiff::derivative(func_rv, autodiff::wrt(r), autodiff::at(plab, r));
                autodiff::derivative(func_iv, autodiff::wrt(plab), autodiff::at(plab, r));
                autodiff::derivative(func_iv, autodiff::wrt(r), autodiff::at(plab, r));
                autodiff::derivative(func_rs, autodiff::wrt(plab), autodiff::at(plab, r));
                autodiff::derivative(func_rs, autodiff::wrt(r), autodiff::at(plab, r));
                autodiff::derivative(func_is, autodiff::wrt(plab), autodiff::at(plab, r));
                autodiff::derivative(func_is, autodiff::wrt(r), autodiff::at(plab, r));
        });
    };
#endif

    BENCHMARK_ADVANCED("Stencil")(Catch::Benchmark::Chronometer meter) {
        double r = 0.15;
        double plab = sqrt(pow(tplab + achilles::Constant::mN, 2) - pow(achilles::Constant::mN, 2));
        auto nucleus = std::make_shared<MockNucleus>();
        REQUIRE_CALL(*nucleus, NNucleons())
            .LR_RETURN((AA))
            .TIMES(AT_LEAST(1));
        achilles::CooperPotential potential;

        meter.measure([&]() {
                potential.derivative_p(nucleus.get(), plab, r);
                potential.derivative_r(nucleus.get(), plab, r);
        });
    };

#endif // CATCH_CONFIG_ENABLE_BENCHMARKING
}

TEST_CASE("CooperPotential::Schroedinger::EDAD1 Values", "[Potential]") {
    constexpr double tplab = 100;
    constexpr size_t AA = 12;

    SECTION("Values are consistent") {
        // Results from the fortran code for tplab = 100, r = 0.15, A = 12
        constexpr double rvector = -27.159583, ivector = -11.257281;

        auto nucleus = std::make_shared<MockNucleus>();
        REQUIRE_CALL(*nucleus, NNucleons())
            .LR_RETURN((AA))
            .TIMES(AT_LEAST(1));

        achilles::SchroedingerPotential potential(5);

        double r = 0.15;
        double plab = sqrt(pow(tplab + achilles::Constant::mN, 2) - pow(achilles::Constant::mN, 2));

        auto vals = potential(nucleus.get(), plab, r);

        // Require to match to 0.01% of Cooper code
        CHECK(vals.rvector == Approx(rvector).epsilon(0.0001));
        CHECK(vals.ivector == Approx(ivector).epsilon(0.0001));
        CHECK(vals.rscalar == Approx(0));
        CHECK(vals.iscalar == Approx(0));
    }
}
