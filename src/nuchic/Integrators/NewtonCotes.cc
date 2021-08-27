#include "nuchic/Integrators/NewtonCotes.hh"

#include "spdlog/spdlog.h"

#include <cmath>
#include <cstdlib>
#include <stdexcept>

using namespace nuchic::Integrators;

NewtonCotes::NewtonCotes(const size_t &order, bool closed) : BaseIntegrator(), m_order(order) {
    if(closed) SetupClosed();
    else SetupOpen();
}

NewtonCotes::NewtonCotes(const size_t &order, bool closed, bool cache) 
    : BaseIntegrator(cache), m_order(order) {
   
    if(closed) SetupClosed();
    else SetupOpen();
}

NewtonCotes::NewtonCotes(const size_t &order, const FunctionD &func, bool closed, bool cache) 
    : BaseIntegrator(func, cache), m_order(order) {

    if(closed) SetupClosed();    
    else SetupOpen();
}

NewtonCotes::NewtonCotes(const size_t &order, const FunctionVD& func, bool closed, bool cache)
    : BaseIntegrator(func, cache), m_order(order) {

    if(closed) SetupClosed();
    else SetupOpen();
}

void NewtonCotes::SetupClosed() {
    if(m_order < 2 || m_order > 4) 
        throw std::runtime_error("Closed Newton-Cotes order must be between 2 and 4.");

    if(m_order == 2) {
        baseFunc = [&](const double &a, const double &h) -> double {
            return h/3.0 * (Function(a) + 4*Function(a+h) + Function(a+2*h));
        };
        errFunc = [&](const double &a, const double &h) -> double {
            return h/2.0 * (Function(a) + Function(a+h)); 
        };
    } else if(m_order == 3) {
        baseFunc = [&](const double &a, const double &h) -> double {
            return 3.0*h/8.0*(Function(a)+3*Function(a+h)+3*Function(a+2*h)+Function(a+3*h));
        };
        errFunc = [&](const double &a, const double &h) -> double {
            return h/3.0*(Function(a)+4*Function(a+h)+Function(a+2*h));
        };
    } else {
        baseFunc = [&](const double &a, const double &h) -> double {
            return 2*h/45.0*(7*Function(a)+32*Function(a+h)+12*Function(a+2*h)
                             +32*Function(a+3*h)+7*Function(a+4*h));
        };
        errFunc = [&](const double &a, const double &h) -> double {
            return 3.0*h/8.0*(Function(a)+3*Function(a+h)+3*Function(a+2*h)+Function(a+3*h));
        };
    }
}

void NewtonCotes::SetupOpen() {
    if(m_order < 3 || m_order > 5) 
        throw std::runtime_error("Open Newton-Cotes order must be between 3 and 5.");

    if(m_order == 3) {
        baseFunc = [&](const double &a, const double &h) -> double {
            return 3.0*h/2.0*(Function(a+h) + Function(a+2*h));
        };
        errFunc = [&](const double &a, const double &h) -> double {
            return 2.0*h*Function(a+h); 
        };
    } else if(m_order == 4) {
        baseFunc = [&](const double &a, const double &h) -> double {
            return 4.0*h/3.0*(2*Function(a+h)-Function(a+2*h)+2*Function(a+3*h));
        };
        errFunc = [&](const double &a, const double &h) -> double {
            return 3.0*h/2.0*(Function(a+h) + Function(a+2*h));
        };
    } else {
        baseFunc = [&](const double &a, const double &h) -> double {
            return 5*h/24.0*(11*Function(a+h)+Function(a+2*h)+Function(a+3*h)+11*Function(a+4*h));
        };
        errFunc = [&](const double &a, const double &h) -> double {
            return 4.0*h/3.0*(2*Function(a+h)-Function(a+2*h)+2*Function(a+3*h));
        };
    }
}

double NewtonCotes::Integrate(const double &a, const double &b, double &err, double&) {
    double h = (b-a)/static_cast<double>(m_order);
    double hErr = (b-a)/static_cast<double>(m_order-1);

    double result = baseFunc(a, h);
    err = std::abs(result - errFunc(a, hErr));

    spdlog::debug("NewtonCotes: I = {} +/- {}", result, err);
    return result;
}

std::vector<double> NewtonCotes::IntegrateVec(const double&, const double&, double&, double&) {
    throw std::runtime_error("Not implemented yet!");
}
