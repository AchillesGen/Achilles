#include "nuchic/Vegas.hh"
#include "nuchic/Utilities.hh"
#include "nuchic/FourVector.hh"

#include "nuchic/Histogram.hh"

#include "yaml-cpp/yaml.h"

#include <vector>
#include <chrono>
#include <iostream>
#include <cmath>

bool histfill = false;
nuchic::Histogram hist(500, 0, 1000, "NumericalDelta");

double delta(double x) {
    constexpr double width2 = 10;
    return 1/sqrt(width2*M_PI)*exp(-x*x/width2);
}

double f(const std::vector<double> &x, const double &wgt) {
    constexpr double mn = 938;
    constexpr double kf = 225;
    constexpr double Ek = 1108;

    double phi = 2*M_PI*x[0];
    double cost = 2*x[1]-1;
    double sint = sqrt(1-cost*cost);
    double Ekp = Ek*x[2];
    double p = kf*x[3];
    double phi2 = 2*M_PI*x[4];
    double cost2 = 2*x[5]-1;
    double sint2 = sqrt(1-cost2*cost2);

    double kp = Ekp;
    double kpx = kp*sint*cos(phi);
    double kpy = kp*sint*sin(phi);
    double kpz = kp*cost;

    double px = p*sint2*cos(phi2);
    double py = p*sint2*sin(phi2);
    double pz = p*cost2;

    double ppx = px-kpx;
    double ppy = py-kpy;
    double ppz = pz+Ek-kpz;
    double Epp = sqrt(ppx*ppx+ppy*ppy+ppz*ppz+mn*mn);

    // Delta function
    // double Ep = 4*kf*x[6];

    // Analytic:
    double Ep = mn+Ek-Ekp-Epp;
    Ep = Ep > 0 ? Ep : 0;
    
    // Brent:
    // auto func = [&](double y){ return mn+Ek-Ekp-Epp-y; };
    // nuchic::Brent brent(func);
    // double Ep = 0;
    // try {
    //     Ep = brent.CalcRoot(0, 1000);
    // } catch (const std::domain_error &e) {
    // }

    double weight = pow(2*M_PI, 2)*4*Ek*Ekp*kp*Ep*kf*p;

    // Delta function
    // weight *= delta(-Ep+mn+Ek-Ekp-Epp)*4*kf;

    if(histfill) hist.Fill(Ep, weight*wgt);

    return weight;
}

int main() {
    YAML::Node node = YAML::Load(R"node(
Vegas:
   seed: 123456789
   iterations: 10
   evaluations: 100000
    )node");

    YAML::Node node2 = YAML::Load(R"node(
Vegas:
   iterations: 20
   evaluations: 10000000
    )node");


    nuchic::AdaptiveMap map(7);
    nuchic::Vegas vegas(map, node["Vegas"]);
    vegas(f);
    vegas.Clear();
    vegas.Set(node2["Vegas"]);
    histfill = true;
    vegas(f);

    hist.Save(hist.GetName());
}
