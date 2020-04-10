#include <cmath>

#include "nuchic/Integrators/Adzint.hh"
#include "spdlog/spdlog.h"

using namespace nuchic::Integrator;

Adzint::Adzint(const int &iacta_, const int &iactb_, bool cache)
    : AdaptiveIntegrator(cache), iacta(iacta_), iactb(iactb_) {}

Adzint::Adzint(const int &iacta_, const int &iactb_, const FunctionS &func, bool cache)
    : AdaptiveIntegrator(func, cache), iacta(iacta_), iactb(iactb_) {}

Adzint::Adzint(const int &iacta_, const int &iactb_, const FunctionV &func, bool cache)
    : AdaptiveIntegrator(func, cache), iacta(iacta_), iactb(iactb_) {}

const double Adzint::sml = 1e-12;
const double Adzint::huge = 1e20;
const double Adzint::tiny = 1e-20;

std::vector<double> Adzint::IntegrateVec(const double&, const double&,
                                         const double&, const double&) {
    throw std::runtime_error("Not implemented!");
}

double Adzint::Integrate(const double &a, const double &b,
                         const double &rerr, const double &aerr) {
    res = 0; ers = 0; ier = 0;
    double ddx = b -a;

    if(ddx == 0) return 0;
    if(ddx < 0) throw std::runtime_error("Integration method(Adzint) only works for b < a");
    if(ddx <= sml) return Function(a+ddx/2)*ddx;

    numint = 3;
    double dx = ddx / static_cast<double>(numint);

    for(size_t i = 0; i < numint; ++i) {
        if(i == 0) {
            u[0] = a;
            if(iacta == 0) fu[0] = Function(a);
            else fa = Function(a+dx/2);
        } else {
            u[i] = v[i-1];
            fu[i] = fv[i-1];
        }
        if(i == numint - 1) {
            v[i] = b;
            if(iactb == 0) Function(b);
            else {
                ib = i;
                fb = Function(b-dx/2);
            }
        } else {
            v[i] = a+dx*static_cast<double>(i+1);
            fv[i] = Function(v[i]);
        }
        Privcal(i, result[i], err[i]);
    }

    for(size_t i = 0; i < numint; ++i) {
        res += result[i];
        ers += err[i];
    }

    const double fac0 = 2.0, facMax=fac0*maxint;
    bool numltmax = true;
    double target = 0;
    while(ier==0 && ers>target && numltmax) {
        size_t numold = numint;
        double facnum = sqrt(static_cast<double>(numold)*fac0);
        do {
            for(size_t i = 0; i < numold; ++i) 
                if(err[i]*facnum > target) Privspl(i);
            if(facnum > facMax) { numltmax=false; break; }
        } while(numold==numint && static_cast<bool>(facnum*=fac0));
        target = std::abs(aerr) + std::abs(rerr*res);
        spdlog::debug("target = {}, aerr = {}, rerr = {}, res = {}", target, aerr, rerr, res);
    }

    if(ier == 1)
        spdlog::warn("Maxint reached (target, error): {}, {}", target, ers);
    spdlog::debug("Result (value, target, error): {}, {}, {}", res, target, ers);

    return res;
}

std::pair<double, double> Adzint::Sglint(const int &iact, const double &f1, const double &f2, 
                                         const double &f3, const double &dx) {
    double tem = dx*(4*f1+3*f2+2*f3)/9;
    double er = dx*(4*f1-6*f2+2*f3)/9;

    if(iact == 2) {
        double t1 = f2 - f1;
        double t2 = f3 - f2;
        if(t1*t2 <= 0) return std::make_pair(tem, er);
        double t12 = t1*t1;
        if(std::abs(f3)*huge < t12) return std::make_pair(tem, er);
        double cc = log(t2/t1)/log(2);
        if(cc <= 1.0) return std::make_pair(tem, er);
        double t3 = t2 - t1;
        double bb = t12/t3;
        double aa = (f1*f3 - f2*f2)/t3;
        double tmp = dx*(aa+bb*pow(4,cc)/(cc+1));
        er = tem - tmp;
        tem = tmp;
    }

    return std::make_pair(tem, er);
}

void Adzint::Privcal(const size_t &i, double &tempResult, double &tempErr) {
    double dx = v[i] - u[i];

    double tem, er;
    if(i == 0 && iacta > 0) {
        fw[i] = fa;
        fa = Function(u[i] + dx/4);
        auto tmp = Sglint(iacta, fa, fw[i], fv[i], dx);
        tem = tmp.first; er = tmp.second;
    } else if(i == ib && iactb > 0) {
        fw[i] = fb;
        fb = Function(v[i] - dx/4);
        auto tmp = Sglint(iactb, fb, fw[i], fu[i], dx);
        tem = tmp.first; er = tmp.second;
    } else {
        double w = (u[i] + v[i]) / 2;
        fw[i] = Function(w);
        tem = dx*(fu[i] + 4*fw[i] + fv[i])/6;
        er = dx*(fu[i] - 2*fw[i] + fv[i])/12;
    }
    tempResult = tem;
    tempErr = std::abs(er);
}

void Adzint::Privspl(const size_t &i) {
    if(numint >= maxint) {
        ier = 1;
        return;
    }

    if(i == ib) ib = numint;
    u[numint] = (u[i]+v[i])/2;
    v[numint] = v[i];

    fu[numint]=fw[i];
    fv[numint]=fv[i];

    v[i] = u[numint];
    fv[i] = fu[numint];

    double oldres = result[i];
    double olderr = err[i];

    Privcal(i, result[i], err[i]);
    Privcal(numint, result[numint], err[numint]);

    double delres = result[i] + result[numint] - oldres;
    res += delres;

    double gooderr = std::abs(delres);

    ers += gooderr - olderr;
    double sumerr = err[i] + err[numint];
    double fac;
    if(sumerr > tiny) fac = gooderr / sumerr;
    else fac = 1.0;

    err[i] *= fac;
    err[numint] *= fac;

    numint++;
}
