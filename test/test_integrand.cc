#include "Achilles/Integrand.hh"
#include "catch2/catch.hpp"

double dummy_func(const std::vector<double> &, double) {
    return 1.0;
}

class DummyMapper : public achilles::Mapper<double> {
  public:
    void GeneratePoint(std::vector<double> &point, const std::vector<double> &rans) override {
        std::copy(rans.begin(), rans.end(), point.begin());
    }
    double GenerateWeight(const std::vector<double> &point, std::vector<double> &rans) override {
        std::copy(point.begin(), point.end(), rans.begin());
        return 1.0;
    }
    size_t NDims() const override { return 1; }
};

TEST_CASE("YAML encoding / decoding Channel", "[multichannel]") {
    achilles::Channel<double> channel;
    achilles::AdaptiveMap map(1, 10);
    channel.integrator = achilles::Vegas(map, achilles::VegasParams{});
    channel.mapping = std::make_unique<DummyMapper>();

    YAML::Node node;
    node["Channel"] = channel;
    auto channel2 = node["Channel"].as<achilles::Channel<double>>();

    CHECK(channel.integrator.Grid().Dims() == channel2.integrator.Grid().Dims());
    CHECK(channel.integrator.Grid().Bins() == channel2.integrator.Grid().Bins());
    CHECK(channel.integrator.Grid().Hist() == channel2.integrator.Grid().Hist());
}

TEST_CASE("YAML encoding / decoding Integrand", "[multichannel]") {
    achilles::Integrand<double> integrand(dummy_func);
    for(size_t i = 0; i < 10; ++i) {
        achilles::Channel<double> channel;
        achilles::AdaptiveMap map(1, 10);
        channel.integrator = achilles::Vegas(map, achilles::VegasParams{});
        channel.mapping = std::make_unique<DummyMapper>();
        integrand.AddChannel(std::move(channel));
    }

    YAML::Node node;
    node["Integrand"] = integrand;
    auto integrand2 = node["Integrand"].as<achilles::Integrand<double>>();

    CHECK(integrand.NChannels() == integrand2.NChannels());
    for(size_t i = 0; i < integrand.NChannels(); ++i) {
        auto grid1 = integrand.GetChannel(i).integrator.Grid();
        auto grid2 = integrand2.GetChannel(i).integrator.Grid();
        CHECK(grid1.Dims() == grid2.Dims());
        CHECK(grid1.Bins() == grid2.Bins());
        CHECK(grid1.Hist() == grid2.Hist());
    }
}
