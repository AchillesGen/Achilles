#include "nuchic/ThreeVector.hh"
#include "nuchic/Rambo.hh"
#include "spdlog/spdlog.h"

nuchic::Rambo::Rambo(size_t nout_) : nout{nout_} {
    double pi2log = log(M_PI/2);
    Z[2] = pi2log;
    for(size_t k = 3; k <= NMAX-nin; ++k) 
        Z[k] = Z[k-1]+pi2log-2*log(static_cast<double>(k-2));
    for(size_t k = 3; k <= NMAX-nin; ++k)
        Z[k] = Z[k]-log(static_cast<double>(k-1));
}

void nuchic::Rambo::GenerateWeight(const std::vector<FourVector> &mom,
                                   const std::vector<double> &mass) {
    FourVector sump(0, 0, 0, 0);
    for(size_t i = 0; i < nin; ++i) sump += mom[i];
    double ET = std::sqrt(sump.M2());
    weight = 1;
    MassiveWeight(mom, mass, ET);
    weight *= exp((2*static_cast<double>(nout)-4)*log(ET)+Z[nout])
              / pow(2*M_PI, static_cast<double>(nout)*3-4);
}

void nuchic::Rambo::GeneratePoint(std::vector<FourVector> &mom,
                                  const std::vector<double> &mass,
                                  const std::vector<double> &rans) {
    FourVector sump(0, 0, 0, 0);
    for(size_t i = 0; i < nin; ++i) sump += mom[i];
    double ET = std::sqrt(sump.M2());

    FourVector R;
    for(size_t i = nin; i < nin+nout; ++i) {
        double ctheta = 2*rans[4*(i-nin)] - 1;
        double stheta = std::sqrt(1-ctheta*ctheta);
        double phi = 2*M_PI*rans[1+4*(i-nin)];
        double Q = -log(rans[2+4*(i-nin)]*rans[3+4*(i-nin)]);
        mom[i] = FourVector(Q*stheta*cos(phi), Q*stheta*sin(phi), Q*ctheta, Q);
        R += mom[i];
    }

    double RMAS = std::sqrt(R.M2());
    ThreeVector B = -R.Vec3()/RMAS;
    double G = R.E()/RMAS;
    double A = 1.0/(1.0+G);
    double X = ET/RMAS;

    for(size_t i = nin; i < nin+nout; ++i) {
        double e = mom[i].E();
        double BQ = B*mom[i].Vec3();
        mom[i] = X*FourVector(mom[i].Vec3()+B*(e+A*BQ), G*e+BQ);
    }

    weight = 1.;
    MassivePoint(mom, mass, ET);
}

void nuchic::Rambo::MassiveWeight(const std::vector<FourVector> &mom,
                                  const std::vector<double> &mass, double ET) {
    static constexpr size_t itmax = 6;
    double accu = ET * std::numeric_limits<double>::epsilon();

    double xmt = 0;
    std::vector<double> E(nin+nout);
    std::vector<double> p2(nin+nout);

    for(size_t i = nin; i < nin+nout; ++i) {
        xmt += mass[i];
        p2[i] = mom[i].E()*mom[i].E();
    }
    double x = 1/std::sqrt(1 - xmt*xmt/ET/ET);

    // Massive particles: Rescale their momenta by a common factor x
    // Loop to calculate x
    for(size_t iter = 0;;) {
        double f0 = -ET;
        double g0 = 0.;
        double x2 = x*x;
        for(size_t i = nin; i < nin+nout; ++i) {
            E[i] = std::sqrt(x2*p2[i]);
            f0 += E[i];
            g0 += p2[i]/E[i];
        }
        if(std::abs(f0) < accu) break;
        if(++iter > itmax) break;
        x -= f0/(x*g0);
    }

    // Calculate momenta + weight
    double wt2 = 1;
    double wt3 = 0;
    for(size_t i = nin; i < nin+nout; ++i) {
        double v = mom[i].P();
        wt2 *= v/mom[i].E();
        wt3 += v*v/mom[i].E();
    }
    x = 1/x;
    weight = exp((2*static_cast<double>(nout)-3)*log(x)+log(wt2/wt3*ET));
}

void nuchic::Rambo::MassivePoint(std::vector<FourVector> &mom,
                                 const std::vector<double> &mass, double ET) const {
    static constexpr size_t itmax = 6;
    double accu = ET * std::numeric_limits<double>::epsilon();

    double xmt = 0;
    std::vector<double> E(nin+nout);
    std::vector<double> p2(nin+nout);

    for(size_t i = nin; i < nin + nout; ++i) {
        xmt += mass[i];
        p2[i] = mom[i].E()*mom[i].E();
    }

    double x = std::sqrt(1 - xmt*xmt/ET/ET);

    // Massive particles: Rescale their momenta by a common factor x
    // Loop to calculate x
    for(size_t iter = 0;;) {
        double f0 = -ET;
        double g0 = 0.;
        double x2 = x*x;
        for(size_t i = nin; i < nin+nout; ++i) {
            E[i] = std::sqrt(mass[i]*mass[i]+x2*p2[i]);
            f0 += E[i];
            g0 += p2[i]/E[i];
        }
        if(std::abs(f0) < accu) break;
        if(++iter > itmax) break;
        x -= f0/(x*g0);
    }

    // Construct Momenta
    for(size_t i = nin; i < nin+nout; ++i) 
        mom[i] = FourVector(x*mom[i].Vec3(), E[i]);
}
