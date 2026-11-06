// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef SYMPLECTIC_INTEGRATORS_HH
#define SYMPLECTIC_INTEGRATORS_HH

#include <functional>
#include <utility>

#include "Achilles/Potential.hh"
#include "Achilles/ThreeVector.hh"

namespace achilles {

/// The extended-phase-space coupling below adds q directly to p, so the two
/// cannot carry a mass dimension here: the integrator works on plain numbers,
/// q in fm and p in MeV. Callers convert at the boundary.
using PSVector = ThreeVectorT<0>;

struct PSState {
    PSVector q, p, x, y;

    PSState() = default;
    PSState(PSVector q_, PSVector p_) : q{q_}, p{p_}, x{q_}, y{p_} {}
};

namespace details {

template <size_t N> struct is_even {
    static constexpr bool value = !(N % 2);
};

} // namespace details

class SymplecticIntegrator {
  public:
    using dHamiltonian = std::function<PSVector(const PSVector &, const PSVector &,
                                                std::shared_ptr<achilles::Potential>)>;
    using PhaseSpace = std::pair<PSVector, PSVector>;
    SymplecticIntegrator() = default;
    SymplecticIntegrator(PSState state, std::shared_ptr<achilles::Potential> pot, dHamiltonian dHdr,
                         dHamiltonian dHdp, double omega)
        : m_omega{std::move(omega)}, m_state{std::move(state)}, m_dHdr{std::move(dHdr)},
          m_dHdp{std::move(dHdp)}, m_pot{std::move(pot)} {}
    SymplecticIntegrator(PSVector q, PSVector p, std::shared_ptr<achilles::Potential> pot,
                         dHamiltonian dHdr, dHamiltonian dHdp, double omega)
        : m_omega{std::move(omega)}, m_dHdr{std::move(dHdr)}, m_dHdp{std::move(dHdp)},
          m_pot{std::move(pot)} {
        m_state = PSState(q, p);
    }

    PSState State() const { return m_state; }
    PSState &State() { return m_state; }
    void Initialize(const PSVector &, const PSVector &);
    PSVector Q() const { return m_state.q; }
    PSVector P() const { return m_state.p; }

    dHamiltonian dHdr() const { return m_dHdr; }
    dHamiltonian &dHdr() { return m_dHdr; }

    dHamiltonian dHdp() const { return m_dHdp; }
    dHamiltonian &dHdp() { return m_dHdp; }

    template <size_t N> void Step(double);

  private:
    void HamiltonianA(double);
    void HamiltonianB(double);
    void Coupling(double);

    double m_omega{1};
    PSState m_state;
    dHamiltonian m_dHdr, m_dHdp;
    std::shared_ptr<Potential> m_pot;
};

template <> inline void SymplecticIntegrator::Step<2>(double time_step) {
    HamiltonianA(time_step / 2);
    HamiltonianB(time_step / 2);
    Coupling(time_step);
    HamiltonianB(time_step / 2);
    HamiltonianA(time_step / 2);
}

template <size_t order> void SymplecticIntegrator::Step(double time_step) {
    static_assert(order % 2 == 0, "SymplecticIntegrator: Order must be an even number");

    // TODO: Make constexpr when c++26 is available
    double gamma = 1.0 / (2 - pow(2, 1.0 / (static_cast<double>(order) + 1.0)));
    Step<order - 2>(gamma * time_step);
    Step<order - 2>((1 - 2 * gamma) * time_step);
    Step<order - 2>(gamma * time_step);
}

} // namespace achilles

#endif
