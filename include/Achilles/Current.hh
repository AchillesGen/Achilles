#pragma once

#include "achilles/ComplexFmt.hh"

#include <array>
#include <complex>
#include <stdexcept>
#include <vector>

namespace achilles {

class VCurrent {
  public:
    VCurrent() = default;
    VCurrent(std::array<std::complex<double>, 4> current) : m_current{current} {}
    VCurrent(std::vector<std::complex<double>> cur) {
        if(cur.size() != 4) throw std::runtime_error("Invalid current size!");
        m_current = std::array<std::complex<double>, 4>{cur[0], cur[1], cur[2], cur[3]};
    }
    std::complex<double> &operator[](size_t lorentz) { return m_current[lorentz]; }
    const std::complex<double> &operator[](size_t lorentz) const { return m_current[lorentz]; }
    std::complex<double> Dot(const VCurrent &other) const {
        return m_current[0] * other.m_current[0] - m_current[1] * other.m_current[1] -
               m_current[2] * other.m_current[2] - m_current[3] * other.m_current[3];
    }
    constexpr size_t size() const { return m_current.size(); }
    constexpr auto begin() { return m_current.begin(); }
    constexpr auto end() { return m_current.end(); }

    VCurrent operator+=(const VCurrent &other) {
        m_current[0] += other.m_current[0];
        m_current[1] += other.m_current[1];
        m_current[2] += other.m_current[2];
        m_current[3] += other.m_current[3];
        return *this;
    }

    bool operator==(const VCurrent &other) const { return m_current == other.m_current; }

  private:
    std::array<std::complex<double>, 4> m_current;
};

inline std::complex<double> operator*(const VCurrent &c1, const VCurrent &c2) {
    return c1.Dot(c2);
}

inline std::vector<VCurrent>
ToCurrentVector(const std::vector<std::vector<std::complex<double>>> &cur) {
    std::vector<VCurrent> result(cur.size());
    for(size_t i = 0; i < cur.size(); ++i) { result[i] = VCurrent(cur[i]); }
    return result;
}

} // namespace achilles
