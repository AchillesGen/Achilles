#include "catch2/catch.hpp"
#include "mock_classes.hh"
#include "nuchic/Potential.hh"
#include "nuchic/Particle.hh"
#include <iostream>

nuchic::PotentialVals<double> stencil5(std::function<nuchic::PotentialVals<double>(double)> f, double x, double h) {
    auto fp2h = f(x + 2*h);
    auto fph = f(x + h);
    auto fm2h = f(x - 2*h);
    auto fmh = f(x - h);

    nuchic::PotentialVals<double> results{};
    double den = 12*h;
    results.rvector = (-fp2h.rvector + 8*fph.rvector - 8*fmh.rvector + fm2h.rvector)/den;
    results.ivector = (-fp2h.ivector + 8*fph.ivector - 8*fmh.ivector + fm2h.ivector)/den;
    results.rscalar = (-fp2h.rscalar + 8*fph.rscalar - 8*fmh.rscalar + fm2h.rscalar)/den;
    results.iscalar = (-fp2h.iscalar + 8*fph.iscalar - 8*fmh.iscalar + fm2h.iscalar)/den;
    return results;
}

TEST_CASE("CooperPotential::EDAD1 Values", "[Potential]") {
    constexpr double tplab = 100;
    constexpr size_t AA = 12;

    SECTION("Values and Derivatives are consistent") {
        autodiff::dual r = 0.15;
    
        // Results from the fortran code for tplab = 100, r = 0.15, A = 12
        constexpr double rvector = 301.65297587389080, ivector = -84.374377568119442;
        constexpr double rscalar = -384.29852972736001, iscalar = 89.69329129030613;

        auto nucleus = std::make_shared<MockNucleus>(); 
        REQUIRE_CALL(*nucleus, NNucleons())
            .LR_RETURN((AA))
            .TIMES(AT_LEAST(1));

        nuchic::CooperPotential potential(nucleus);

        autodiff::dual plab = sqrt(pow(tplab + nuchic::Constant::mN, 2) - pow(nuchic::Constant::mN, 2));
        nuchic::PotentialVals<autodiff::dual> vals = potential.evaluate<autodiff::dual>(plab, r);

        auto func_rv = [&](const autodiff::dual &x, const autodiff::dual &y){ return potential.evaluate<autodiff::dual>(x, y).rvector; };
        auto func_iv = [&](const autodiff::dual &x, const autodiff::dual &y){ return potential.evaluate<autodiff::dual>(x, y).ivector; };
        auto func_rs = [&](const autodiff::dual &x, const autodiff::dual &y){ return potential.evaluate<autodiff::dual>(x, y).rscalar; };
        auto func_is = [&](const autodiff::dual &x, const autodiff::dual &y){ return potential.evaluate<autodiff::dual>(x, y).iscalar; };
        double drvdp = autodiff::derivative(func_rv, autodiff::wrt(plab), autodiff::at(plab, r));
        double drvdr = autodiff::derivative(func_rv, autodiff::wrt(r), autodiff::at(plab, r));
        double divdp = autodiff::derivative(func_iv, autodiff::wrt(plab), autodiff::at(plab, r));
        double divdr = autodiff::derivative(func_iv, autodiff::wrt(r), autodiff::at(plab, r));
        double drsdp = autodiff::derivative(func_rs, autodiff::wrt(plab), autodiff::at(plab, r));
        double drsdr = autodiff::derivative(func_rs, autodiff::wrt(r), autodiff::at(plab, r));
        double disdp = autodiff::derivative(func_is, autodiff::wrt(plab), autodiff::at(plab, r));
        double disdr = autodiff::derivative(func_is, autodiff::wrt(r), autodiff::at(plab, r));

        // auto [drvdp, drvdr] = autodiff::derivatives(vals.rvector, autodiff::wrt(plab, r));
        // auto [divdp, divdr] = autodiff::derivatives(vals.ivector, autodiff::wrt(plab, r));
        // auto [drsdp, drsdr] = autodiff::derivatives(vals.rscalar, autodiff::wrt(plab, r));
        // auto [disdp, disdr] = autodiff::derivatives(vals.iscalar, autodiff::wrt(plab, r));

        CHECK(vals.rvector == Approx(rvector));
        CHECK(vals.ivector == Approx(ivector));
        CHECK(vals.rscalar == Approx(rscalar));
        CHECK(vals.iscalar == Approx(iscalar));

        std::cout << drvdp << " " << drvdr << "\n";
        std::cout << divdp << " " << divdr << "\n";
        std::cout << drsdp << " " << drsdr << "\n";
        std::cout << disdp << " " << disdr << "\n";

        std::function<nuchic::PotentialVals<double>(double)> func_p = [&](double x){ return potential.evaluate<double>(x, static_cast<double>(r));};
        std::function<nuchic::PotentialVals<double>(double)> func_r = [&](double x){ return potential.evaluate<double>(static_cast<double>(plab), x);};
        auto stencilp = stencil5(func_p, static_cast<double>(plab), 0.01);
        auto stencilr = stencil5(func_r, static_cast<double>(r), 0.01);

        CHECK(drvdp == Approx(stencilp.rvector));
        CHECK(drvdr == Approx(stencilr.rvector));
        CHECK(divdp == Approx(stencilp.ivector));
        CHECK(divdr == Approx(stencilr.ivector));
        CHECK(drsdp == Approx(stencilp.rscalar));
        CHECK(drsdr == Approx(stencilr.rscalar));
        CHECK(disdp == Approx(stencilp.iscalar));
        CHECK(disdr == Approx(stencilr.iscalar));

    }

#if defined(CATCH_CONFIG_ENABLE_BENCHMARKING)

    BENCHMARK_ADVANCED("Autodiff")(Catch::Benchmark::Chronometer meter) {
        autodiff::dual r = 0.15;
        auto nucleus = std::make_shared<MockNucleus>(); 
        REQUIRE_CALL(*nucleus, NNucleons())
            .LR_RETURN((AA))
            .TIMES(AT_LEAST(1));
        nuchic::CooperPotential potential(nucleus);
        autodiff::dual plab = sqrt(pow(tplab + nuchic::Constant::mN, 2) - pow(nuchic::Constant::mN, 2));

        auto func_rv = [&](const autodiff::dual &x, const autodiff::dual &y){ return potential.evaluate<autodiff::dual>(x, y).rvector; };
        auto func_iv = [&](const autodiff::dual &x, const autodiff::dual &y){ return potential.evaluate<autodiff::dual>(x, y).ivector; };
        auto func_rs = [&](const autodiff::dual &x, const autodiff::dual &y){ return potential.evaluate<autodiff::dual>(x, y).rscalar; };
        auto func_is = [&](const autodiff::dual &x, const autodiff::dual &y){ return potential.evaluate<autodiff::dual>(x, y).iscalar; };

        meter.measure([&]() {
                potential.evaluate<autodiff::dual>(plab, r);
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

    BENCHMARK_ADVANCED("Stencil")(Catch::Benchmark::Chronometer meter) {
        double r = 0.15;
        double plab = sqrt(pow(tplab + nuchic::Constant::mN, 2) - pow(nuchic::Constant::mN, 2));
        auto nucleus = std::make_shared<MockNucleus>(); 
        REQUIRE_CALL(*nucleus, NNucleons())
            .LR_RETURN((AA))
            .TIMES(AT_LEAST(1));
        nuchic::CooperPotential potential(nucleus);
        std::function<nuchic::PotentialVals<double>(double)> func_p = [&](double x){ return potential.evaluate<double>(x, r); };
        std::function<nuchic::PotentialVals<double>(double)> func_r = [&](double x){ return potential.evaluate<double>(plab, x); };

        meter.measure([&]() {
                stencil5(func_p, plab, 0.01);
                stencil5(func_r, r, 0.01);
        });
    };

#endif // CATCH_CONFIG_ENABLE_BENCHMARKING
}

