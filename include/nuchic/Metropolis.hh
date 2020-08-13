#ifndef METROPOLIS_HH
#define METROPOLIS_HH

#include <array>
#include <limits>

namespace nuchic {

class Particle;

class Metropolis {
    public:
        Metropolis() = default;
        bool IsInelastic(const double &ran, const double &energy, bool isSame) const;

    private:
        static constexpr std::array<double, 8> energies = {335, 410, 510, 660, 840, 1160, 1780, 3900};
        static constexpr std::array<double, 8> fiiInel = {0.07, 0.20, 0.31, 0.43, 0.58, 0.65, 0.69, 0.69};
        static constexpr std::array<double, 8> fijInel = {0.04, 0.07, 0.15, 0.27, 0.37, 0.36, 0.35, 0.35};
        static constexpr std::array<double, 8> fpi = {1.0, 1.0, 1.0, 1.0, 0.97, 0.80, 0.44, 0.44};
        static constexpr std::array<double, 8> Aii = {0.1, 0.9, 2.7, 9.0, 14.3, 19.2, std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
        static constexpr std::array<double, 8> Bii = {0, 0, 0, 0, 0, 0, 0, 0};
        static constexpr std::array<double, 8> Aij = {2.2, 1.8, 2.3, 8.8, 15.0, 29.4, std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
        static constexpr std::array<double, 8> Bij = {-1.0, -1.1, -0.7, -0.2, 0, 0, 0, 0};
};

//TODO: come up with better names for these
enum class PionInteractionMode {
    same, different, zero
};

class PionInteractions {
    public:
        PionInteractions() = default;
        double xsec(const Particle&, const Particle&) const;
        bool IsInelastic(const double &ran, const Particle&, const Particle&) const;

    private:
        PionInteractionMode GetMode(const Particle&, const Particle&) const;

        static constexpr std::array<double, 8> energies = {49, 85, 128, 184, 250, 350, 540, 1300};
        static constexpr std::array<double, 8> sigmaii = {16, 50, 114, 200, 110, 51, 20, 30};
        static constexpr std::array<double, 8> sigmaij = {15, 21, 43, 66, 44, 23, 22, 30};
        static constexpr std::array<double, 8> sigmaabs = {20, 32, 45, 36, 18, 0, 0, 0};
        static constexpr std::array<double, 8> fiiInel = {0, 0, 0, 0.03, 0.06, 0.16, 0.30, 0.88};
        static constexpr std::array<double, 8> fiiCEX = {0, 0, 0, 0, 0, 0, 0, 0};
        static constexpr std::array<double, 8> fijInel = {0.45, 0.57, 0.62, 0.64, 0.62, 0.56, 0.58, 0.94};
        static constexpr std::array<double, 8> fijCEX = {1.0, 1.0, 1.0, 0.95, 0.89, 0.72, 0.51, 0.06};
        static constexpr std::array<double, 8> f0Inel = {0.42, 0.36, 0.36, 0.37, 0.40, 0.50, 0.59, 0.94};
        static constexpr std::array<double, 8> f0CEX = {1.0, 1.0, 1.0, 0.90, 0.84, 0.67, 0.50, 0.05};
        static constexpr std::array<double, 8> fpi = {1.0, 1.0, 1.0, 1.0, 1.0, 0.98, 0.91, 0.24};
        static constexpr std::array<double, 8> Aii = {3.2, 2.2, 1.9, 2.2, 2.6, 3.0, 3.0, 3.0};
        static constexpr std::array<double, 8> Bii = {-1.8, -2.1, -1.5, -0.3, 2.0, 4.0, 4.0, 4.0};
        static constexpr std::array<double, 8> Aij = {1.1, 1.9, 2.2, 2.2, 2.0, 2.7, 3.0, 3.0};
        static constexpr std::array<double, 8> Bij = {0.8, 0.7, 0.8, 1.0, 1.4, 2.6, 3.6, 4.0};
        static constexpr std::array<double, 8> A0 = {3.4, 2.1, 1.9, 2.1, 2.5, 3.0, 3.0, 3.0};
        static constexpr std::array<double, 8> B0 = {-1.8, -2.0, -1.4, 0, 1.7, 4.0, 4.0, 4.0};
};

}

#endif
