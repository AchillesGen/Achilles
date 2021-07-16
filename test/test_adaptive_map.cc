#include "catch2/catch.hpp" 

#include "nuchic/AdaptiveMap.hh"

#include "catch_utils.hh"
#include <iostream>
#include <sstream>

TEST_CASE("Adaptive Map Construction", "[vegas]") {
    nuchic::AdaptiveMap2 map(2, 4);

    CHECK(map.Dims() == 2);
    CHECK(map.Bins() == 4);

    const std::vector<double> edges{0.0, 0.25, 0.5, 0.75, 1.0};
    CHECK(edges == map.Edges(0));
    CHECK(edges == map.Edges(1));
}

TEST_CASE("Numbers generated according to Adaptive Map", "[vegas]") {
    constexpr size_t ndims = 2;
    constexpr size_t nbins = 10;
    nuchic::AdaptiveMap2 map(ndims, nbins);

    SECTION("Uniform Distribution"){
        auto rans = GENERATE(take(100, randomVector(ndims)));
        const auto output = rans;
        auto wgt = map(rans);

        // Check that the jacobian is 1 for a uniform distribution
        CHECK(wgt == Approx(1.0));

        // Check that the input number is the same as the output numbers;
        CHECK_THAT(rans, Catch::Matchers::Approx(rans));
    }

    SECTION("Non-Uniform Distribution") {
        auto rans = GENERATE(take(100, randomVector(ndims)));
        const auto input = rans;
        const std::vector<double> newEdges{0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0,
                                           0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0};
        const std::vector<double> jacobian{2,   1,   1,   1,   1,   1,   1,   1, 0.5, 0.5};
        map.Hist() = newEdges;
        auto wgt = map(rans);

        // Check that the jacobian is correct
        double jac = 1.0;
        for(size_t i = 0; i < ndims; ++i) {
            const auto position = input[i] * nbins;
            const auto idx = static_cast<size_t>(position);
            jac *= jacobian[idx];
            const auto val = newEdges[idx] + (position - static_cast<double>(idx)) * map.width(i, idx);
            CHECK(rans[i] == Approx(val));
        }

        CHECK(wgt == Approx(jac));
    }
}

TEST_CASE("Adaptive Map Histogram Updates", "[vegas]") {
    SECTION("Adapting the map") {
        constexpr size_t ndims = 2;
        constexpr size_t nbins = 10;
        nuchic::AdaptiveMap2 map(ndims, nbins);
        const auto data = GENERATE(take(1, randomVector(ndims*nbins, 0, 100)));

        // No change if alpha = 0
        std::array<std::vector<double>, ndims> edges;
        for(size_t i = 0; i < ndims; ++i) edges[i] = map.Edges(i);
        map.Adapt(0, data); 
        for(size_t i = 0; i < ndims; ++i)
            CHECK_THAT(map.Edges(i), Catch::Matchers::Approx(edges[i]));

        // Change if alpha != 0
        map.Adapt(1.5, data); 
        for(size_t i = 0; i < ndims; ++i) {
            CHECK(std::is_sorted(map.Hist().begin()+static_cast<int>(i*(nbins+1)),
                                 map.Hist().begin()+static_cast<int>((i+1)*(nbins+1))));
            for(size_t j = 1; j < nbins; ++j) {
                CHECK(map.Edges(i)[j] != Approx(edges[i][j]));
                CHECK((map.Edges(i)[j] >= 0 && map.Edges(i)[j] <= 1));
            }
        }
    }

    SECTION("Splitting the map") {
        nuchic::AdaptiveMap2 map(2, 2);

        map.Split();
        CHECK(map.Dims() == 2);
        CHECK(map.Bins() == 4);
        const std::vector<double> edges{0.0, 0.25, 0.5, 0.75, 1.0};
        CHECK(edges == map.Edges(0));
        CHECK(edges == map.Edges(1));

        map = nuchic::AdaptiveMap2(2, 2);
        map.Split(nuchic::AdaptiveMapSplit::third);
        CHECK(map.Dims() == 2);
        CHECK(map.Bins() == 6);
        const std::vector<double> edges2{0.0, 0.5/3.0, 1.0/3.0, 1.5/3.0, 2.0/3.0, 2.5/3.0, 1.0};
        CHECK_THAT(map.Edges(0), Catch::Matchers::Approx(edges2));
        CHECK_THAT(map.Edges(1), Catch::Matchers::Approx(edges2));

        map = nuchic::AdaptiveMap2(2, 2);
        map.Split(nuchic::AdaptiveMapSplit::quarter);
        CHECK(map.Dims() == 2);
        CHECK(map.Bins() == 8);
        const std::vector<double> edges3{0.0, 0.5/4, 1.0/4, 1.5/4, 2.0/4, 2.5/4, 3.0/4, 3.5/4, 1.0};
        CHECK_THAT(map.Edges(0), Catch::Matchers::Approx(edges3));
        CHECK_THAT(map.Edges(1), Catch::Matchers::Approx(edges3));
    }
}

TEST_CASE("Serializing / Deserializing Adaptive Map", "[vegas]") {
    nuchic::AdaptiveMap2 map(4, 4);

    std::stringstream data;
    map.Serialize(data);

    nuchic::AdaptiveMap2 map2;
    map2.Deserialize(data);

    CHECK(map.Dims() == map2.Dims());
    CHECK(map.Bins() == map2.Bins());

    for(size_t i = 0; i < map.Dims(); ++i)
        CHECK_THAT(map.Edges(i), Catch::Matchers::Approx(map2.Edges(i)));
}
