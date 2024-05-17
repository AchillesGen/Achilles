#ifndef RANDOM_HH
#define RANDOM_HH

#include <iterator>
#include <memory>

#include "Achilles/Randutils.hh"

namespace achilles {

class Random {
  public:
    static Random Instance() {
        static Random rand;
        return rand;
    }

    void Seed(unsigned int seed) { m_rng->seed(seed); }

    void Generate(std::vector<double> &vec) {
        m_rng->generate<std::uniform_real_distribution>(vec);
    }

    template <typename T, size_t N> void Generate(std::array<T, N> &array, T low = 0, T high = 1) {
        m_rng->generate(array, low, high);
    }

    template <typename T> T Uniform(T low, T high) { return m_rng->uniform(low, high); }

    template <typename T> T Pick(const std::vector<T> &vec) { return m_rng->pick(vec); }

    template <typename T> std::size_t SelectIndex(const T &array) {
        return m_rng->variate<std::size_t, std::discrete_distribution>(array);
    }

    template <typename T> void Shuffle(std::vector<T> &vec) {
        return m_rng->shuffle(vec.begin(), vec.end());
    }

    std::size_t SelectIndex(const std::vector<double> &array) {
        return m_rng->variate<std::size_t, std::discrete_distribution>(array.begin(), array.end());
    }

    template <typename T>
    void Sample(size_t n, const std::vector<T> &array, std::vector<T> &result) {
        std::sample(array.begin(), array.end(), std::back_inserter(result), n, m_rng->engine());
    }

    void SaveState(std::ostream &os) { os << m_rng->engine(); }
    void LoadState(std::istream &is) { is >> m_rng->engine(); }

  private:
    Random() { m_rng = std::make_shared<randutils::mt19937_rng>(); }
    std::shared_ptr<randutils::mt19937_rng> m_rng;
};

} // namespace achilles

#endif
