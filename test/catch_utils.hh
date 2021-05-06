#ifndef CATCH_UTILS_HH
#define CATCH_UTILS_HH

// Inspired from Catch2/examples/300-Gen-OwnGenerator.cpp

#include <string>
#include <random>
#include "catch2/catch.hpp"

#include "nuchic/FourVector.hh"

class RandomNucleusGenerator : public Catch::Generators::IGenerator<std::string> {
    static constexpr auto& chars = "abcdefghijklmnopqrstuvwxyz"
                                   "ABCDEFGHIJKLMNOPQRSTUVWXYZ";

    size_t m_max_nucleons;
    std::mt19937 m_rand{std::random_device{}()};
    std::uniform_int_distribution<std::string::size_type> m_dist;
    std::uniform_int_distribution<size_t> m_dist_nucleons;
    std::string current_string;

    public:
        RandomNucleusGenerator(size_t max_nucleons):
            m_max_nucleons{max_nucleons},
            m_rand(std::random_device{}()),
            m_dist(0, sizeof(chars) - 2),
            m_dist_nucleons(1, m_max_nucleons) {

            static_cast<void>(next());
        };

        std::string const& get() const override;
        bool next() override {
            current_string = "";
            current_string += std::to_string(m_dist_nucleons(m_rand));
            for(size_t i = 0; i < 2; ++i)
                current_string += chars[m_dist(m_rand)];
            return true;
        }
};

class RandomMomentumGenerator : public Catch::Generators::IGenerator<nuchic::FourVector> {
    double m_maxVal;
    std::mt19937 m_rand{std::random_device{}()};
    std::uniform_real_distribution<double> m_dist;
    nuchic::FourVector current_momentum{};

    public:
        RandomMomentumGenerator(double maxVal):
            m_maxVal{maxVal},
            m_rand(std::random_device{}()),
            m_dist(0, m_maxVal) {

            static_cast<void>(next());
        };

        nuchic::FourVector const& get() const override;
        bool next() override {
            double px = m_dist(m_rand);
            double py = m_dist(m_rand);
            double pz = m_dist(m_rand);
            double e = m_dist(m_rand);
            current_momentum = nuchic::FourVector(e, px, py, pz);
            return true;
        }
};

#endif
