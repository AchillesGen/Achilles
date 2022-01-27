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
    double m_mass;
    std::mt19937 m_rand{std::random_device{}()};
    std::uniform_real_distribution<double> m_dist;
    nuchic::FourVector current_momentum{};

    public:
        RandomMomentumGenerator(double maxVal, double mass = 0):
            m_maxVal{maxVal}, m_mass(mass),
            m_rand(std::random_device{}()),
            m_dist(0, m_maxVal) {

            static_cast<void>(next());
        };

        nuchic::FourVector const& get() const override;
        bool next() override {
            double pmag = m_dist(m_rand);
            double cost = 2*m_dist(m_rand)/m_maxVal - 1;
            double sint = sqrt(1-cost*cost);
            double phi = 2*M_PI*m_dist(m_rand)/m_maxVal;
            double e = sqrt(pmag*pmag+m_mass*m_mass);
            double px = pmag*sint*cos(phi);
            double py = pmag*sint*sin(phi);
            double pz = pmag*cost;
            current_momentum = nuchic::FourVector(e, px, py, pz);
            return true;
        }
};

inline nuchic::FourVector const& RandomMomentumGenerator::get() const {
    return current_momentum;
}
// This helper function provides a nicer UX when instantiating the generator
// Notice that it returns an instance of GeneratorWrapper<std::string>, which
// is a value-wrapper around std::unique_ptr<IGenerator<std::string>>.
inline Catch::Generators::GeneratorWrapper<nuchic::FourVector> randomMomentum(double max, double mass = 0) {
    return Catch::Generators::GeneratorWrapper<nuchic::FourVector>(std::unique_ptr<Catch::Generators::IGenerator<nuchic::FourVector>>(new RandomMomentumGenerator(max, mass)));
}

class RandomVectorGenerator : public Catch::Generators::IGenerator<std::vector<double>> {
    double m_minVal, m_maxVal;
    size_t m_size;

    std::mt19937 m_rand{std::random_device{}()};
    std::uniform_real_distribution<double> m_dist;
    std::vector<double> current_vector{};

    public:
        RandomVectorGenerator(size_t size, double minVal=0.0, double maxVal=1.0):
            m_minVal{minVal}, m_maxVal{maxVal}, m_size{size},
            m_rand(std::random_device{}()),
            m_dist(m_minVal, m_maxVal) {

            current_vector.resize(m_size);
            static_cast<void>(next());
        };

        std::vector<double> const& get() const override;
        bool next() override {
            for(size_t i = 0; i < m_size; ++i) current_vector[i] = m_dist(m_rand);

            return true;
        }
};

inline std::vector<double> const & RandomVectorGenerator::get() const {
    return current_vector;
}

inline Catch::Generators::GeneratorWrapper<std::vector<double>> randomVector(size_t size, double min=0, double max=1) {
    return Catch::Generators::GeneratorWrapper<std::vector<double>>(std::unique_ptr<Catch::Generators::IGenerator<std::vector<double>>>(new RandomVectorGenerator(size, min, max)));
}


#endif
