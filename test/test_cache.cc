#include "Achilles/MultiChannel.hh"
#include "Achilles/Process.hh"
#include "Achilles/Vegas.hh"
#include "catch2/catch.hpp"
#include "catch_utils.hh"
#include "mock_classes.hh"

#include "Achilles/Random.hh"
#include "Achilles/Unweighter.hh"

#include <sstream>

TEST_CASE("Save/Load Random State", "[Cache]") {
    auto instance = achilles::Random::Instance();
    instance.Seed(42);
    std::stringstream ss;
    instance.SaveState(ss);

    std::vector<double> vec(10);
    instance.Generate(vec);

    instance.LoadState(ss);
    std::vector<double> vec2(10);
    instance.Generate(vec2);

    CHECK(vec == vec2);
}

TEST_CASE("Save/Load Statistics", "[Cache]") {
    SECTION("Percentile") {
        achilles::Percentile percentile(0.9);
        auto data = GENERATE(take(2, randomVector(100)));
        for(const auto &d : data) { percentile.Add(d); }
        std::stringstream ss;
        percentile.SaveState(ss);

        achilles::Percentile percentile2(0.9);
        percentile2.LoadState(ss);
        CHECK(percentile.Get() == percentile2.Get());
    }

    SECTION("StatsData") {
        achilles::StatsData stats;
        auto data = GENERATE(take(2, randomVector(100)));
        for(const auto &d : data) { stats += d; }
        std::stringstream ss;
        stats.SaveState(ss);

        achilles::StatsData stats2;
        stats2.LoadState(ss);
        CHECK(stats == stats2);
    }
}

TEMPLATE_TEST_CASE("Save/Load Unweighters", "[Cache]", achilles::NoUnweighter,
                   achilles::PercentileUnweighter) {
    TestType unweighter(YAML::Load("percentile: 90"));
    auto data = GENERATE(take(2, randomVector(100)));
    for(const auto &d : data) { unweighter.AddEvent(d); }

    std::stringstream ss;
    unweighter.SaveState(ss);

    TestType unweighter2(YAML::Load("percentile: 90"));
    unweighter2.LoadState(ss);
    CHECK(unweighter.MaxValue() == unweighter2.MaxValue());

    auto val = GENERATE(take(1, random(0.0, 1.0)));

    // Need to ensure that the state is the same since a random number is used to accept or not
    std::stringstream ss2;
    achilles::Random::Instance().SaveState(ss2);
    auto accept1 = unweighter.AcceptEvent(val);
    achilles::Random::Instance().LoadState(ss2);
    auto accept2 = unweighter2.AcceptEvent(val);
    CHECK(accept1 == accept2);
}

TEST_CASE("Save/Load Vegas Parameters", "[Cache]") {
    achilles::VegasParams params{1000, 1, 1e-8, 1e-8, 2, 100};

    std::stringstream ss;
    achilles::SaveState(ss, params);

    achilles::VegasParams params2;
    achilles::LoadState(ss, params2);
    CHECK(params == params2);
}

TEST_CASE("Save/Load Vegas Summary", "[Cache]") {
    achilles::VegasSummary summary;
    for(size_t i = 0; i < 3; ++i) {
        achilles::StatsData stats;
        auto data = GENERATE(take(2, randomVector(100)));
        for(const auto &d : data) { stats += d; }
        summary.results.push_back(stats);
        summary.sum_results += stats;
    }

    std::stringstream ss;
    achilles::SaveState(ss, summary);

    achilles::VegasSummary summary2;
    achilles::LoadState(ss, summary2);
    CHECK(summary == summary2);
}

TEST_CASE("Save/Load Adaptive Map", "[Cache]") {
    achilles::AdaptiveMap map(2, 10);
    auto data = GENERATE(take(2, randomVector(100)));
    map.Adapt(1.0, data);

    std::stringstream ss;
    map.SaveState(ss);

    achilles::AdaptiveMap map2(2, 10);
    map2.LoadState(ss);
    CHECK(map.Hist() == map2.Hist());
}

double test_func_cache(const std::vector<double> &x, double wgt) {
    return 3.0 / 2.0 * (x[0] * x[0] + x[1] * x[1]) * wgt;
}

TEST_CASE("Save/Load Vegas", "[Cache]") {
    achilles::Vegas vegas({2, 10}, {1000, 1, 1e-2, 1e-2, 2, 2});
    vegas.SetVerbosity(0);
    vegas.Optimize(test_func_cache);
    std::stringstream ss;
    vegas.SaveState(ss);

    achilles::Vegas vegas2;
    vegas2.LoadState(ss);

    CHECK(vegas == vegas2);
}

TEST_CASE("Save/Load MultiChannel Parameters", "[Cache]") {
    achilles::MultiChannelParams params{1000, 10, 1e-2, 1, 0.25, 1e-5, 0};

    std::stringstream ss;
    achilles::SaveState(ss, params);

    achilles::MultiChannelParams params2;
    achilles::LoadState(ss, params2);
    CHECK(params == params2);
}

TEST_CASE("Save/Load MultiChannel Summary", "[Cache]") {
    achilles::MultiChannelSummary summary;
    summary.best_weights = {0.3, 0.4, 0.1, 0.2};
    for(size_t i = 0; i < 3; ++i) {
        achilles::StatsData stats;
        auto data = GENERATE(take(2, randomVector(100)));
        for(const auto &d : data) { stats += d; }
        summary.results.push_back(stats);
        summary.sum_results += stats;
    }

    std::stringstream ss;
    achilles::SaveState(ss, summary);

    achilles::MultiChannelSummary summary2;
    achilles::LoadState(ss, summary2);
    CHECK(summary == summary2);
}

constexpr double s0 = -10.0;
constexpr double s1 = 10.0;

double test_func_exp_cache(const std::vector<double> &x, double wgt) {
    double sms0 = x[0] - s0;
    double sms1 = x[0] - s1;
    return (2.0 * std::exp(-sms0 * sms0) + std::exp(-sms1 * sms1)) / std::sqrt(std::acos(-1.0)) /
           3.0 * wgt;
}

class DoubleMapperCache : public achilles::Mapper<double> {
  public:
    DoubleMapperCache(size_t channel) : m_channel{std::move(channel)} {}
    void GeneratePoint(std::vector<double> &point, const std::vector<double> &rans) override {
        double s = std::tan(std::acos(-1.0) * (rans[0] - 0.5));
        if(m_channel == 0)
            s += s0;
        else
            s += s1;
        point[0] = s;
    }
    double GenerateWeight(const std::vector<double> &point, std::vector<double> &rans) override {
        const double sms0 = point[0] - s0;
        const double sms1 = point[0] - s1;
        rans.resize(1);
        if(m_channel == 1) {
            rans[0] = std::atan(sms1) / std::acos(-1.0) + 0.5;
            return 1.0 / (1.0 + sms1 * sms1) / std::acos(-1.0);
        }
        rans[0] = std::atan(sms0) / std::acos(-1.0) + 0.5;
        return 1.0 / (1.0 + sms0 * sms0) / std::acos(-1.0);
    }
    size_t NDims() const override { return 1; }
    YAML::Node ToYAML() const override { return YAML::Node(); }

  private:
    size_t m_channel;
};

TEST_CASE("Save/Load MultiChannel", "[Cache]") {
    achilles::Integrand<double> integrand(test_func_exp_cache);
    for(size_t i = 0; i < 2; ++i) {
        achilles::Channel<double> channel;
        channel.mapping = std::make_unique<DoubleMapperCache>(i);
        achilles::AdaptiveMap map(channel.mapping->NDims(), 50);
        channel.integrator = achilles::Vegas(map, achilles::VegasParams{});
        integrand.AddChannel(std::move(channel));
    }
    achilles::MultiChannel multi(1, integrand.NChannels(), {1000, 10, 1e-2, 1, 0.25, 1e-5, 0});
    multi.Optimize(integrand);
    std::stringstream ss;
    multi.SaveState(ss);

    achilles::MultiChannel multi2;
    multi2.LoadState(ss);

    CHECK(multi == multi2);
}

TEST_CASE("Save/Load ProcessInfo", "[Cache]") {
    achilles::ProcessInfo info;
    info.m_leptonic = {achilles::PID::electron(), {achilles::PID::electron()}};
    info.m_hadronic = {{achilles::PID::proton()}, {achilles::PID::proton()}};
    info.m_spectator = {};

    std::stringstream ss;
    info.SaveState(ss);

    achilles::ProcessInfo info2;
    info2.LoadState(ss);

    CHECK(info == info2);
}

TEST_CASE("Save/Load Process", "[Cache]") {
    auto unweighter = std::make_unique<MockUnweighter>();
    REQUIRE_CALL(*unweighter, AddEvent(10)).TIMES(1000);
    REQUIRE_CALL(*unweighter, MaxValue()).TIMES(1).RETURN(10);
    REQUIRE_CALL(*unweighter, SaveState(trompeloeil::_)).TIMES(2);
    achilles::ProcessInfo info;
    info.m_leptonic = {achilles::PID::electron(), {achilles::PID::electron()}};
    info.m_hadronic = {{achilles::PID::proton()}, {achilles::PID::proton()}};
    info.m_spectator = {};

    achilles::Process process(info, std::move(unweighter));
    for(size_t i = 0; i < 1000; ++i) process.AddWeight(10);

    std::stringstream ss;
    process.SaveState(ss);
    process.SaveState(ss);
    std::cout << ss.str() << std::endl;

    auto unweighter2 = std::make_unique<MockUnweighter>();
    REQUIRE_CALL(*unweighter2, LoadState(trompeloeil::_)).TIMES(1);
    REQUIRE_CALL(*unweighter2, MaxValue()).TIMES(1).RETURN(10);
    achilles::Process process2(info, std::move(unweighter2));
    process2.LoadState(ss);

    CHECK(process.Info() == process2.Info());
    CHECK(process.MaxWeight() == process2.MaxWeight());
    CHECK(process.TotalCrossSection() == process2.TotalCrossSection());
}
