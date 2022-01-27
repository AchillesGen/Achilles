#ifndef RANDOM_HH
#define RANDOM_HH

#include <memory>

#include "nuchic/Randutils.hh"

namespace nuchic {

class Random {
    public:
        static Random Instance() {
            static Random rand;
            return rand;
        }

        void Seed(unsigned int seed) {
            m_rng -> seed(seed);
        }

        void Generate(std::vector<double>& vec) {
            m_rng -> generate<std::uniform_real_distribution>(vec);
        }

        template<typename T, size_t N>
        void Generate(std::array<T, N> &array, T low=0, T high=1) {
            m_rng -> generate(array, low, high);
        }

        template<typename T>
        T Uniform(T low, T high) {
            return m_rng -> uniform(low, high);
        }

        template<typename T>
        T Pick(const std::vector<T> &vec) {
            return m_rng -> pick(vec);
        }

        template<typename T>
        std::size_t SelectIndex(const T &array) {
            return m_rng -> variate<std::size_t, std::discrete_distribution>(array);
        }

        std::size_t SelectIndex(const std::vector<double> &array) {
            return m_rng -> variate<std::size_t, std::discrete_distribution>(array.begin(), array.end());
        }

    private:
        Random() {
            m_rng = std::make_shared<randutils::mt19937_rng>();
        }
        std::shared_ptr<randutils::mt19937_rng> m_rng;
};

}

#endif
