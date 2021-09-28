#include "nuchic/Potential.hh"
#include "nuchic/Nucleus.hh"
#include <iostream>

using nuchic::PotentialVals;
using nuchic::Potential;
using nuchic::CooperPotential;
using nuchic::SchroedingerPotential;

constexpr std::array<double, 8*8> CooperPotential::pt1;
constexpr std::array<double, 8*8> CooperPotential::pt2;
constexpr std::array<double, 8*6> CooperPotential::pt3;

PotentialVals Potential::stencil5(std::function<nuchic::PotentialVals(double)> f, double x, double h) const {
    auto fp2h = f(x + 2*h);
    auto fph = f(x + h);
    auto fm2h = f(x - 2*h);
    auto fmh = f(x - h);

    nuchic::PotentialVals results{};
    double den = 12*h;
    results.rvector = (-fp2h.rvector + 8*fph.rvector - 8*fmh.rvector + fm2h.rvector)/den;
    results.ivector = (-fp2h.ivector + 8*fph.ivector - 8*fmh.ivector + fm2h.ivector)/den;
    results.rscalar = (-fp2h.rscalar + 8*fph.rscalar - 8*fmh.rscalar + fm2h.rscalar)/den;
    results.iscalar = (-fp2h.iscalar + 8*fph.iscalar - 8*fmh.iscalar + fm2h.iscalar)/den;
    return results;
}

PotentialVals Potential::stencil5second(std::function<nuchic::PotentialVals(double)> f,
                                        double x, double h) const {
    auto fp2h = f(x + 2*h);
    auto fph = f(x + h);
    auto fm2h = f(x - 2*h);
    auto fmh = f(x - h);
    auto f0 = f(x);

    nuchic::PotentialVals results{};
    double den = 12*pow(h,2);
    results.rvector = (-fp2h.rvector + 16*fph.rvector - 30*f0.rvector +16*fmh.rvector - fm2h.rvector)/den;
    results.ivector = (-fp2h.ivector + 16*fph.ivector - 30*f0.ivector +16*fmh.ivector - fm2h.ivector)/den;
    results.rscalar = (-fp2h.rscalar + 16*fph.rscalar - 30*f0.rscalar +16*fmh.rscalar - fm2h.rscalar)/den;
    results.iscalar = (-fp2h.iscalar + 16*fph.iscalar - 30*f0.iscalar +16*fmh.iscalar - fm2h.iscalar)/den;
    return results;
}

std::array<PotentialVals, 3> Potential::stencil5all(std::function<nuchic::PotentialVals(double)> f,
                                                    double x, double h) const {
    auto fp2h = f(x + 2*h);
    auto fph = f(x + h);
    auto fm2h = f(x - 2*h);
    auto fmh = f(x - h);
    auto f0 = f(x);

    std::array<nuchic::PotentialVals, 3> results{};
    results[0] = f0;

    // First derivative
    double den = 12*h;
    results[1].rvector = (-fp2h.rvector + 8*fph.rvector - 8*fmh.rvector + fm2h.rvector)/den;
    results[1].ivector = (-fp2h.ivector + 8*fph.ivector - 8*fmh.ivector + fm2h.ivector)/den;
    results[1].rscalar = (-fp2h.rscalar + 8*fph.rscalar - 8*fmh.rscalar + fm2h.rscalar)/den;
    results[1].iscalar = (-fp2h.iscalar + 8*fph.iscalar - 8*fmh.iscalar + fm2h.iscalar)/den;

    // Second derivative
    den = 12*pow(h,2);
    results[2].rvector = (-fp2h.rvector + 16*fph.rvector - 30*f0.rvector +16*fmh.rvector - fm2h.rvector)/den;
    results[2].ivector = (-fp2h.ivector + 16*fph.ivector - 30*f0.ivector +16*fmh.ivector - fm2h.ivector)/den;
    results[2].rscalar = (-fp2h.rscalar + 16*fph.rscalar - 30*f0.rscalar +16*fmh.rscalar - fm2h.rscalar)/den;
    results[2].iscalar = (-fp2h.iscalar + 16*fph.iscalar - 30*f0.iscalar +16*fmh.iscalar - fm2h.iscalar)/den;

    return results;
}

std::unique_ptr<Potential> nuchic::SquareWellPotential::Construct(std::shared_ptr<Nucleus>& nuc,
                                                                  const YAML::Node&) {
    return std::make_unique<SquareWellPotential>(nuc); 
}

PotentialVals nuchic::SquareWellPotential::operator()(const double&, const double &r) const {
    PotentialVals result;
    constexpr double potential_shift = 8;
    result.rvector = -(sqrt(Constant::mN*Constant::mN + pow(m_nucleus->FermiMomentum(r), 2))
        - Constant::mN + potential_shift);

    return result;
}

std::unique_ptr<Potential> nuchic::WiringaPotential::Construct(std::shared_ptr<Nucleus>& nuc,
                                                               const YAML::Node &node) {
    double r0 = node["r0"].as<double>();
    return std::make_unique<WiringaPotential>(nuc, r0);
}

PotentialVals nuchic::WiringaPotential::operator()(const double &plab, const double &radius) const {
    const double rho = m_nucleus -> Rho(radius);
    const double rho_ratio = rho/m_rho0;
    const double alpha = 15.52*rho_ratio + 24.93*pow(rho_ratio, 2);
    const double beta = -116*rho_ratio;
    const double lambda = (3.29 - 0.373*rho_ratio)*nuchic::Constant::HBARC;

    PotentialVals results{};
    results.rvector = alpha + beta/(1+pow(plab/lambda, 2));
    return results;
}

std::unique_ptr<Potential> CooperPotential::Construct(std::shared_ptr<Nucleus>& nuc,
                                                      const YAML::Node&) {
    return std::make_unique<CooperPotential>(nuc);
}

PotentialVals CooperPotential::evaluate(const double &plab, const double &radius) const {
    const auto tplab = sqrt(plab*plab + pow(nuchic::Constant::mN, 2)) - nuchic::Constant::mN;
    const auto aa = static_cast<double>(m_nucleus -> NNucleons());
    const auto wt = aa * nuchic::Constant::AMU;
    const auto ee = tplab;
    const auto acb = cbrt(aa);
    const auto caa = aa / (aa + 20);
    const auto y = caa, y2 = y*y, y3 = y*y2, y4 = y2*y2;
    const auto el = ee+wp;
    const auto wp2 = wp*wp;
    const auto wt2 = wt*wt;
    const auto pcm = sqrt(wt2*(el*el-wp2)/(wp2+wt2+2.0*wt*el));
    const auto epcm = sqrt(wp2+pcm*pcm);
    const auto etcm = sqrt(wt2+pcm*pcm);
    const auto sr = epcm + etcm;
    const auto e = 1000.0/epcm;
    const auto x = e;
    const auto x2 = x*x;
    const auto x3 = x*x2;
    const auto x4 = x2*x2;
    const auto recv = (etcm / sr);
    const auto recs = (wt / sr);
    constexpr double cv1b = 1.0, cv2b = 1.0, cs1b = 1.0, cs2b = 1.0;
    constexpr double av1b = 0.7, av2b = 0.7, as1b = 0.7, as2b = 0.7;
    const auto sumr = -100. * (Data(1, 1)+Data(1, 2)*x+Data(1, 3)*x2+Data(1, 4)*x3+Data(1, 5)*x4
                                         +Data(1, 6)*y+Data(1, 7)*y2+Data(1, 8)*y3+Data(20,1)*y4);
    const auto rv1 =   cv1b * (Data(2, 1)+Data(2, 2)*x+Data(2, 3)*x2+Data(2, 4)*x3+Data(2, 5)*x4
                                         +Data(2, 6)*y+Data(2, 7)*y2+Data(2, 8)*y3+Data(21,1)*y4);
    const auto av1 =   av1b * (Data(3, 1)+Data(3, 2)*x+Data(3, 3)*x2+Data(3, 4)*x3+Data(3, 5)*x4
                                         +Data(3, 6)*y+Data(3, 7)*y2+Data(3, 8)*y3);
    const auto sumi = -15.0 * (Data(4, 1)+Data(4, 2)*x+Data(4, 3)*x2+Data(4, 4)*x3+Data(4, 5)*x4
                                         +Data(4, 6)*y+Data(4, 7)*y2+Data(4, 8)*y3+Data(20,2)*y4);
    const auto rv2 =   cv2b * (Data(5, 1)+Data(5, 2)*x+Data(5, 3)*x2+Data(5, 4)*x3+Data(5, 5)*x4
                                         +Data(5, 6)*y+Data(5, 7)*y2+Data(5, 8)*y3+Data(21,2)*y4);
    const auto av2 =   av2b * (Data(6, 1)+Data(6, 2)*x+Data(6, 3)*x2+Data(6, 4)*x3+Data(6, 5)*x4
                                         +Data(6, 6)*y+Data(6, 7)*y2+Data(6, 8)*y3);
    const auto diffr = 700. * (Data(7, 1)+Data(7, 2)*x+Data(7, 3)*x2+Data(7, 4)*x3+Data(7, 5)*x4
                                         +Data(7, 6)*y+Data(7, 7)*y2+Data(7, 8)*y3+Data(20,3)*y4);
    const auto rs1 =   cs1b * (Data(8, 1)+Data(8, 2)*x+Data(8, 3)*x2+Data(8, 4)*x3+Data(8, 5)*x4
                                         +Data(8, 6)*y+Data(8, 7)*y2+Data(8, 8)*y3+Data(21,3)*y4);
    const auto as1 =   as1b * (Data(9, 1)+Data(9, 2)*x+Data(9, 3)*x2+Data(9, 4)*x3+Data(9, 5)*x4
                                         +Data(9, 6)*y+Data(9, 7)*y2+Data(9, 8)*y3);
    const auto diffi = -150 * (Data(10, 1)+Data(10, 2)*x+Data(10, 3)*x2+Data(10, 4)*x3+Data(10, 5)*x4
                                          +Data(10, 6)*y+Data(10, 7)*y2+Data(10, 8)*y3+Data(20,4)*y4);
    const auto rs2 =   cs2b * (Data(11, 1)+Data(11, 2)*x+Data(11, 3)*x2+Data(11, 4)*x3+Data(11, 5)*x4
                                          +Data(11, 6)*y+Data(11, 7)*y2+Data(11, 8)*y3+Data(21,4)*y4);
    const auto as2 =   as2b * (Data(12, 1)+Data(12, 2)*x+Data(12, 3)*x2+Data(12, 4)*x3+Data(12, 5)*x4
                                          +Data(12, 6)*y+Data(12, 7)*y2+Data(12, 8)*y3);

    const auto vv = 0.5*(sumr+diffr);
    const auto vs = 0.5*(sumr-diffr);
    const auto wv = 0.5*(sumi+diffi);
    const auto ws = 0.5*(sumi-diffi);

    const auto wv2 = -100.0*(Data(13, 1)+Data(13, 2)*x+Data(13, 3)*x2+Data(13, 4)*x3+Data(13, 5)*x4
                                        +Data(13, 6)*y+Data(13, 7)*y2+Data(13, 8)*y3+Data(20, 5)*y4);
    // const auto rv22 =       (Data(14, 1)+Data(14, 2)*x+Data(14, 3)*x2+Data(14, 4)*x3+Data(14, 5)*x4);
    // const auto av22 =   0.7*(Data(15, 1)+Data(15, 2)*x+Data(15, 3)*x2+Data(15, 4)*x3+Data(15, 5)*x4);

    const auto ws2 =  100.0*(Data(16, 1)+Data(16, 2)*x+Data(16, 3)*x2+Data(16, 4)*x3+Data(16, 5)*x4
                                        +Data(16, 6)*y+Data(16, 7)*y2+Data(16, 8)*y3+Data(20, 6)*y4);
    // const auto rs22 =       (Data(17, 1)+Data(17, 2)*x+Data(17, 3)*x2+Data(17, 4)*x3+Data(17, 5)*x4);
    // const auto as22 =   0.7*(Data(18, 1)+Data(18, 2)*x+Data(18, 3)*x2+Data(18, 4)*x3+Data(18, 5)*x4);

    const auto rv22=rv2;
    const auto av22=av2;
    const auto rs22=rs2;
    const auto as22=as2;

    const auto rva1 = CalcTerm(recv*vv, rv1, acb, av1, radius);
    const auto rva2 = CalcTerm(recv*wv, rv2, acb, av2, radius) + CalcTermSurf(recv*wv2, rv22, acb, av22, radius);
    const auto rsa1 = CalcTerm(recs*vs, rs1, acb, as1, radius);
    const auto rsa2 = CalcTerm(recs*ws, rs2, acb, as2, radius) + CalcTermSurf(recs*ws2, rs22, acb, as22, radius);

    return {rva1, rsa1, rva2, rsa2};
}

std::unique_ptr<Potential> SchroedingerPotential::Construct(std::shared_ptr<Nucleus>& nuc,
                                                            const YAML::Node &node) {
    size_t mode = node["Mode"].as<size_t>();
    return std::make_unique<SchroedingerPotential>(nuc, mode);
}

nuchic::PotentialVals SchroedingerPotential::operator()(const double &plab, const double &radius) const {
    auto potential = stencil5all([&](double r){ return evaluate(plab, r); }, radius, 0.01);
    const double u1 = potential[0].rscalar;
    const double w1 = potential[0].iscalar;
    const double u2 = potential[0].rvector;
    const double w2 = potential[0].ivector;
    const double ud1 = potential[1].rscalar;
    const double udd1 = potential[2].rscalar;
    const double wd1 = potential[1].iscalar;
    const double wdd1 = potential[2].iscalar;
    const double ud2 = potential[1].rvector;
    const double udd2 = potential[2].rvector;
    const double wd2 = potential[1].ivector;
    const double wdd2 = potential[2].ivector;

    const auto tplab = sqrt(plab*plab + pow(nuchic::Constant::mN, 2)) - nuchic::Constant::mN;
    const auto aa = static_cast<double>(m_nucleus -> NNucleons());
    const auto wt = aa * nuchic::Constant::AMU;
    const auto ee = tplab;
    const auto el = ee+wp;
    const auto wp2 = wp*wp;
    const auto wt2 = wt*wt;
    const auto pcm = sqrt(wt2*(el*el-wp2)/(wp2+wt2+2.0*wt*el));
    const auto epcm = sqrt(wp2+pcm*pcm);
    const double etcm = sqrt(wt2+pcm*pcm);
    const double sr = epcm+etcm;

    double redu{}; //this corresponds to opt 5 needs more investigation
    switch(m_mode) {
        case 2:
            redu = epcm*etcm/sr;
            break;
        case 3:
            redu = wp*wt/(wp+wt);
            break;
        case 4:
            redu = epcm;
            break;
        case 5:
            redu = wp;
            break;
    }
    const auto ac = 0; //can be changed if Coulomb corrections are to be included
    const auto dc = 0; //can be changed if Coulomb corrections are to be included
    const auto ddc = 0; //can be changed if Coulomb corrections are to be included
    const auto couf1 = 0;
    const auto couf2 = 0;

    const auto ucr1 = 2*epcm*u2 + 2*wp*u1-pow(u2,2)+pow(w2,2)+pow(u1,2)-pow(w1,2)-2*ac*u2;
    const auto ucrw = 0.5/redu*ucr1;
    const auto uci1 = 2*(epcm*w2+wp*w1-u2*w2+u1*w1-ac*w2);
    const auto uciw = 0.5/redu*uci1;

    // central terms
    const auto aab = epcm+wp+u1-u2-ac;
    const auto adr = ud1-ud2 -dc;
    const auto adi = wd1-wd2;
    const auto addr = udd1-udd2-ddc;
    const auto addi = wdd1-wdd2;
    const auto aai = w1-w2;

    const auto udr1 = -1.0/radius*(adr*aab+adi*aai)/(pow(aab,2)+pow(aai,2));
    const auto udr2 = -0.5*(addr*aab+addi*aai)/(pow(aab,2)+pow(aai,2));
    const auto udr3 = 0.75*((pow(adr,2)-pow(adi,2))
                             *(pow(aab,2)-pow(aai,2))+4.0*adr*adi*aab*aai)
                            /(pow(pow(aab,2)-pow(aai,2),2)+4.0*pow(aab,2)*pow(aai,2));
    const auto udrw = 0.5*hc2/redu*(udr1+udr2+udr3);
    const auto udi1 = -1.0/radius*(adi*aab-adr*aai)/(pow(aab,2)+pow(aai,2));
    const auto udi2 = -0.5*(addi*aab-addr*aai)/(pow(aab,2)+pow(aai,2));
    const auto udi3 = 1.5*(adr*adi*(pow(aab,2)-pow(aai,2))-(pow(adr,2)-pow(adi,2))*aab*aai)
                           /(pow(pow(aab,2)-pow(aai,2),2)+4.0*pow(aab,2)*pow(aai,2));
    const auto udiw = 0.5*hc2/redu*(udi1+udi2+udi3);

    // CALCULATION OF CENTRAL POTENTIAL (WITH DARWIN)
    const auto uer = ucrw+udrw-couf2*0.5*pow(ac,2)/redu+couf1*(epcm/redu)*ac;
    const auto uei = uciw+udiw;

    PotentialVals results{};
    results.rvector=uer;
    results.ivector=uei;
    return results;
}
