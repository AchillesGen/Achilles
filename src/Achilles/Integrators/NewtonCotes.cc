#include "Achilles/Integrators/NewtonCotes.hh"
#include <stdexcept>

using namespace achilles::Integrator;

NewtonCotes::NewtonCotes(const size_t &order, bool closed)
    : QuadratureIntegrator(), m_order(order) {
    if(closed)
        SetupClosed();
    else
        SetupOpen();
}

NewtonCotes::NewtonCotes(const size_t &order, bool closed, bool cache)
    : QuadratureIntegrator(cache), m_order(order) {
    if(closed)
        SetupClosed();
    else
        SetupOpen();
}

NewtonCotes::NewtonCotes(const size_t &order, const FunctionS &func, bool closed, bool cache)
    : QuadratureIntegrator(func, cache), m_order(order) {
    if(closed)
        SetupClosed();
    else
        SetupOpen();
}

NewtonCotes::NewtonCotes(const size_t &order, const FunctionV &func, bool closed, bool cache)
    : QuadratureIntegrator(func, cache), m_order(order) {
    if(closed)
        SetupClosed();
    else
        SetupOpen();
}

void NewtonCotes::SetupClosed() {
    if(m_order < 2 || m_order > 4)
        throw std::runtime_error("Closed Newton-Cotes order must be between 2 and 4.");

    if(m_order == 2) {
        baseFunc = [&](const double &a, const double &h) -> double {
            return h / 3 * (Function(a) + 4 * Function(a + h) + Function(a + 2 * h));
        };
        errFunc = [&](const double &a, const double &h) -> double {
            return h / 2 * (Function(a) + Function(a + h));
        };
    } else if(m_order == 3) {
        baseFunc = [&](const double &a, const double &h) -> double {
            return 3 * h / 8 *
                   (Function(a) + 3 * Function(a + h) + 3 * Function(a + 2 * h) +
                    Function(a + 3 * h));
        };
        errFunc = [&](const double &a, const double &h) -> double {
            return h / 3 * (Function(a) + 4 * Function(a + h) + Function(a + 2 * h));
        };
    } else {
        baseFunc = [&](const double &a, const double &h) -> double {
            return 2 * h / 45 *
                   (7 * Function(a) + 32 * Function(a + h) + 12 * Function(a + 2 * h) +
                    32 * Function(a + 3 * h) + 7 * Function(a + 4 * h));
        };
        errFunc = [&](const double &a, const double &h) -> double {
            return 3 * h / 8 *
                   (Function(a) + 3 * Function(a + h) + 3 * Function(a + 2 * h) +
                    Function(a + 3 * h));
        };
    }
}

void NewtonCotes::SetupOpen() {
    if(m_order < 3 || m_order > 5)
        throw std::runtime_error("Open Newton-Cotes order must be between 3 and 5.");

    if(m_order == 3) {
        baseFunc = [&](const double &a, const double &h) -> double {
            return 3 * h / 2 * (Function(a + h) + Function(a + 2 * h));
        };
        errFunc = [&](const double &a, const double &h) -> double {
            return 2 * h * Function(a + h);
        };
    } else if(m_order == 4) {
        baseFunc = [&](const double &a, const double &h) -> double {
            return 4 * h / 3 *
                   (2 * Function(a + h) - Function(a + 2 * h) + 2 * Function(a + 3 * h));
        };
        errFunc = [&](const double &a, const double &h) -> double {
            return 3 * h / 2 * (Function(a + h) + Function(a + 2 * h));
        };
    } else {
        baseFunc = [&](const double &a, const double &h) -> double {
            return 5 * h / 24 *
                   (11 * Function(a + h) + Function(a + 2 * h) + Function(a + 3 * h) +
                    11 * Function(a + 4 * h));
        };
        errFunc = [&](const double &a, const double &h) -> double {
            return 4 * h / 3 *
                   (2 * Function(a + h) - Function(a + 2 * h) + 2 * Function(a + 3 * h));
        };
    }
}

double NewtonCotes::Integrate(const double &a, const double &b, double &err) {
    double h = (b - a) / static_cast<double>(m_order);
    double hErr = (b - a) / static_cast<double>(m_order - 1);

    double result = baseFunc(a, h);
    err = std::abs(result - errFunc(a, hErr));

    return result;
}

std::vector<double> NewtonCotes::IntegrateVec(const double &, const double &, double &) {
    throw std::runtime_error("Not implemented yet!");
}
