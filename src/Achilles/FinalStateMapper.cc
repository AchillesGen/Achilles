#include "Achilles/FinalStateMapper.hh"
#include "Achilles/Constants.hh"
#include "Achilles/FourVector.hh"
#include "Achilles/Poincare.hh"
#include "Achilles/ThreeVector.hh"
#include "Achilles/Units.hh"
#include "Achilles/Utilities.hh"

#ifdef ACHILLES_SHERPA_INTERFACE
using achilles::SherpaMapper;
using ATOOLS::Vec4D;
#endif
using achilles::FourVector;
using achilles::ThreeBodyMapper;
using achilles::TwoBodyMapper;

void TwoBodyMapper::GeneratePoint(std::vector<FourVector> &mom, const std::vector<double> &rans) {
    // The momentum are given in the following order:
    // 1. Momentum of the initial hadron
    // 2. Momentum of the initial lepton
    // 3. Momentum of all outgoing parts of the leptonic tensor
    // 4. Momentum of all outgoing hadrons
    auto p01 = (mom[0] + mom[1]);
    auto s = p01.M2();
    auto sqrts = sqrt(s);
    auto boostVec = p01.BoostVector();
    auto mom0 = mom[0].Boost(-boostVec);
    Poincare zax(mom0, FourVector(1., 0., 0., 1.));
    auto cosT = dCos * rans[0] - 1;
    auto sinT = sqrt(1 - cosT * cosT);
    auto phi = dPhi * rans[1];
    auto E1 = sqrts / 2 * (1 + s2 / s - s3 / s);
    auto E2 = sqrts / 2 * (1 + s3 / s - s2 / s);
    auto lambda = sqrt(pow(s - s2 - s3, 2) - 4 * s2 * s3);
    auto pCM = lambda / (2 * sqrts);

    mom[2] = {E1, pCM * sinT * cos(phi), pCM * sinT * sin(phi), pCM * cosT};
    mom[3] = {E2, -pCM * sinT * cos(phi), -pCM * sinT * sin(phi), -pCM * cosT};

    zax.RotateBack(mom[2]);
    zax.RotateBack(mom[3]);

    mom[2] = mom[2].Boost(boostVec);
    mom[3] = mom[3].Boost(boostVec);

    // Mapper<achilles::FourVector>::Print(__PRETTY_FUNCTION__, mom, rans);
    // spdlog::trace("  MassCheck: {}", CheckMasses({mom[2], mom[3]}, {s2, s3}));
    // spdlog::trace("  s = {}, lambda = {}", s, lambda);
}

double TwoBodyMapper::GenerateWeight(const std::vector<FourVector> &mom,
                                     std::vector<double> &rans) {
    auto boostVec = (mom[0] + mom[1]).BoostVector();
    auto mom0 = mom[0].Boost(-boostVec);
    auto rotMat = mom0.AlignZ();
    auto p2 = mom[2].Boost(-boostVec).Rotate(rotMat);
    rans[0] = (p2.CosTheta() + 1) / dCos;
    rans[1] = p2.Phi() / dPhi;

    auto pcm = p2.P();
    auto ecm = (mom[0] + mom[1]).M();

    auto factor = pcm / ecm / (16 * M_PI * M_PI);
    auto wgt = 1.0 / dCos / dPhi / factor;
    // Mapper<achilles::FourVector>::Print(__PRETTY_FUNCTION__, mom, rans);
    // spdlog::trace("  ct: {}", p2.CosTheta());
    // spdlog::trace("  pcm: {}", pcm);
    // spdlog::trace("  ecm: {}", ecm);
    // spdlog::trace("  Weight: {}", wgt);

    return wgt;
}

void ThreeBodyMapper::GeneratePoint(std::vector<FourVector> &mom, const std::vector<double> &rans) {
    // Mass and width of "resonance"
    // double res_mass = 1232.;
    // double res_width = 120.;

    auto p01 = (mom[0] + mom[1]);
    auto s = p01.M2();
    auto sqrts = sqrt(s);

    // Minimum and maximum invariant mass of hadrons
    double s23_max = pow(sqrts - sqrt(s4), 2);
    double s23_min = std::max(pow(sqrt(s2) + sqrt(s3), 2), 1E-8);

    // old
    // auto s23 = MassivePropagatorMomenta(res_mass, res_width, s23_min, s23_max, rans[0]);

    // new
    auto s23 = s23_min + (s23_max - s23_min) * rans[0];

    FourVector p23;
    // Should tmass = res_mass?
    TChannelMomenta(mom[0], mom[1], p23, mom[4], s23, s4, 0.0, m_alpha, m_ctmax, m_ctmin, m_amct, 0,
                    rans[1], rans[2]);
    Isotropic2Momenta(p23, s2, s3, mom[2], mom[3], rans[3], rans[4], m_ctmin, m_ctmax);

    // Want to know gammaN invariant mass
    // auto pGamma = mom[1] - mom[4];
    // auto pGammaN = pGamma + mom[0];
    // auto GammaN_mass = pGammaN.M();
    // spdlog::debug("GammaN invariant mass = {}", GammaN_mass);

    // Mapper<achilles::FourVector>::Print(__PRETTY_FUNCTION__, mom, rans);
}

double ThreeBodyMapper::GenerateWeight(const std::vector<FourVector> &mom,
                                       std::vector<double> &rans) {
    // The momentum are given in the following order:
    // 1. Momentum of the initial hadron
    // 2. Momentum of the initial lepton
    // 3. Momentum of all outgoing hadrons
    // 4. Momentum of all outgoing parts of the leptonic tensor

    // Mass and width of "resonance"
    // double res_mass = 1232.;
    // double res_width = 120.;
    auto wt = 1.0;

    auto p01 = (mom[0] + mom[1]);
    auto s = p01.M2();
    auto sqrts = sqrt(s);

    double s23_max = pow(sqrts - sqrt(s4), 2);
    double s23_min = std::max(pow(sqrt(s2) + sqrt(s3), 2), 1E-8);

    auto s23 = (mom[2] + mom[3]).M2();
    rans[0] = (s23 - s23_min) / (s23_max - s23_min);
    wt *= 1. / (s23_max - s23_min);

    // wt *=
    // MassivePropWeight(res_mass,res_width,1,s23_min,s23_max,fabs((mom[2]+mom[3]).M2()),rans[0]);

    wt *= TChannelWeight(mom[0], mom[1], mom[2] + mom[3], mom[4], 0.0, m_alpha, m_ctmax, m_ctmin,
                         m_amct, 0, rans[1], rans[2]);
    wt *= Isotropic2Weight(mom[2], mom[3], rans[3], rans[4], m_ctmin, m_ctmax);

    if(wt != 0.) wt = 1.0 / wt / pow(2. * M_PI, (3 * 3.) - 4.);
    // Mapper<achilles::FourVector>::Print(__PRETTY_FUNCTION__, mom, rans);
    return 1 / wt;
}

double ThreeBodyMapper::MassivePropagatorMomenta(double mass, double width, double smin,
                                                 double smax, double ran) {
    double mass2 = pow(mass, 2);
    auto mw = mass * width;
    auto ymax = atan((smin - mass2) / mw);
    auto ymin = atan((smax - mass2) / mw);
    auto s = mass2 + mw * tan(ymin + ran * (ymax - ymin));

    if(!(s > 0) && !(s < 0) && s != 0) { spdlog::debug("MassivePropMomenta produced a nan"); }

    return s;
}

void ThreeBodyMapper::Isotropic2Momenta(FourVector p, double s1_, double s2_, FourVector &p1,
                                        FourVector &p2, double ran1, double ran2, double ctmin,
                                        double ctmax) {
    auto s = p.M2();
    auto rs = sqrt(fabs(s));
    FourVector p1h;
    p1h.SetE((s + s1_ - s2_) / rs / 2.);
    auto p1m = rs * SqLam(s, s1_, s2_) / 2.;
    auto ct = ctmin + (ctmax - ctmin) * ran1;
    auto st = sqrt(1. - pow(ct, 2));
    auto phi = dPhi * ran2;

    p1h = {p1h.E(), p1m * st * sin(phi), p1m * st * cos(phi), p1m * ct};
    Boost(0, p, p1h, p1);
    p2 = p + (-1.) * p1;

    spdlog::debug("  MassCheck: {}", CheckMasses({p1, p2}, {s1_, s2_}));
}

int ThreeBodyMapper::TChannelMomenta(FourVector p1in, FourVector p2in, FourVector &p1out,
                                     FourVector &p2out, double s1out, double s2out, double t_mass,
                                     double ctexp, double ctmax, double ctmin, double aminct,
                                     int aminctflag, double ran1, double ran2) {
    auto t_mass2 = t_mass * t_mass;
    FourVector pin = p1in + p2in;
    auto s = pin.M2();
    auto sabs = sqrt(fabs(s));
    auto s1in = p1in.M2();
    auto s2in = p2in.M2();

    FourVector p1inh, p1outh;
    p1inh.SetE((s + s1in - s2in) / 2. / sabs);
    auto p1inmass = sabs * SqLam(s, s1in, s2in) / 2.;
    p1inh = {p1inh.E(), 0., 0., p1inmass};

    p1outh.SetE((s + s1out - s2out) / 2. / sabs);
    auto p1outmass = sabs * SqLam(s, s1out, s2out) / 2.;

    auto a = (t_mass2 - s1in - s1out + 2. * p1outh.E() * p1inh.E()) / (2. * p1inmass * p1outmass);
    if(a <= 1.0 + 1.0e-6) a = 1.0 + 1.0e-6;
    if(a < aminct) a = aminct;

    if(fabs(a - ctmax) < 1.e-14) a = ctmax;

    aminct = Tj1(ctexp, a - ctmin, a - ctmax, ran1);
    auto ct = a - aminct;
    double st;
    if(aminctflag == 1)
        st = sqrt(aminct * (1. + ct));
    else
        st = sqrt(1. - pow(ct, 2));
    auto phi = dPhi * ran2;
    p1outh = {p1outh.E(), p1outmass * st * cos(phi), p1outmass * st * sin(phi), p1outmass * ct};

    FourVector help;
    Boost(1, pin, p1in, help);

    Poincare Rot(p1inh, help);
    help = p1outh;
    Rot.Rotate(help);
    Boost(0, pin, help, p1out);

    p2out = pin + (-1.) * p1out;
    spdlog::debug("MassCheck: {}", CheckMasses({p1out, p2out}, {s1out, s2out}));
    return 0;
}

double ThreeBodyMapper::SqLam(double s, double s1_, double s2_) {
    double arg1 = pow(s - s1_ - s2_, 2) - 4. * s1_ * s2_;
    if(arg1 > 0.)
        return sqrt(arg1) / s;
    else
        return 0.;
}

double ThreeBodyMapper::Tj1(double cn, double amcxm, double amcxp, double ran) {
    double ce = 1. - cn;
    double res = 0.;
    if(ce > 1.e-12)
        res = pow(ran * pow(amcxm, ce) + (1. - ran) * pow(amcxp, ce), 1. / ce);
    else {
        if(amcxp > 0.)
            res = exp(ran * log(amcxm) + (1. - ran) * log(amcxp));
        else
            res = -exp(ran * log(-amcxm) + (1. - ran) * log(-amcxp));
    }
    return res;
}

double ThreeBodyMapper::Hj1(double cn, double amcxm, double amcxp) {
    double ce = 1. - cn;
    if(ce > 1.e-12) return (pow(amcxp, ce) - pow(amcxm, ce)) / ce;
    return log(amcxp / amcxm);
}

double ThreeBodyMapper::MassivePropWeight(double mass, double width, int lim, double smin,
                                          double smax, double s, double &ran) {
    double mass2 = mass * mass;
    double mw = mass * width;
    if(lim == 0)
        return mw / (M_PI * ((s - mass2) * (s - mass2) + mw * mw));
    else {
        if((s < smin) || (s > smax) || smin == smax) {
            ran = -1.;
            return 0.;
        }

        double ymax = atan((smin - mass2) / mw);
        double ymin = atan((smax - mass2) / mw);
        double y = atan((s - mass2) / mw);
        ran = (y - ymin) / (ymax - ymin);

        double wt = mw / ((s - mass2) * (s - mass2) + mw * mw);
        wt /= ymin - ymax;

        return wt;
    }
}

double ThreeBodyMapper::TChannelWeight(const FourVector &p1in, const FourVector &p2in,
                                       const FourVector &p1out, const FourVector &p2out,
                                       double t_mass, double ctexp, double ctmax, double ctmin,
                                       double aminct, int, double &ran1, double &ran2) {
    double t_mass2 = t_mass * t_mass;
    FourVector pin = p1in + p2in;
    double s = pin.M2();
    double sabs = sqrt(fabs(s));
    double s1in = p1in.M2();
    double s2in = p2in.M2();
    double s1out = p1out.M2();
    double s2out = p2out.M2();
    if(s1out < 1.e-8) s1out = 0.;
    if(s2out < 1.e-8) s2out = 0.;
    FourVector p1inh, p1outh;
    p1inh.SetE((s + s1in - s2in) / 2. / sabs);
    double p1inmass = sabs * SqLam(s, s1in, s2in) / 2.;
    p1inh = {p1inh.E(), 0., 0., p1inmass};
    p1outh.SetE((s + s1out - s2out) / 2. / sabs);
    double p1outmass = sabs * SqLam(s, s1out, s2out) / 2.;

    double a = (t_mass2 - s1in - s1out + 2. * p1outh.E() * p1inh.E()) / (2. * p1inmass * p1outmass);
    if(a <= 1.0 + 1.0e-6) a = 1.0 + 1.0e-6;
    if(a < aminct) a = aminct;

    FourVector help = p1out;
    Boost(1, pin, help, p1outh);
    help = p1in;
    Boost(1, pin, help, p1inh);

    Poincare Rot(FourVector(1., 0., 0., 1.), p1inh);
    Rot.RotateBack(p1outh);

    double pa1;
    if(fabs(a - ctmax) < 1.e-14)
        pa1 = 0.;
    else
        pa1 = pow(a - ctmax, 1. - ctexp);
    double ct = p1outh.Z() / p1outh.P();
    if(ct < ctmin || ct > ctmax) {
        ran1 = ran2 = -1.;
        return 0.;
    }
    ran1 = (pow(a - ct, 1. - ctexp) - pa1);
    ran1 /= (pow(a - ctmin, 1. - ctexp) - pa1);
    ran2 = asin(p1outh.Y() / p1outh.Pt()) / (dPhi);
    if(p1outh.X() < 0.) ran2 = .5 - ran2;
    if(ran2 < 0.) ran2 += 1.;

    aminct = a - ct;
    double wt =
        2. * sabs / (-pow(aminct, ctexp) * Hj1(ctexp, a - ctmin, a - ctmax) * p1outmass * M_PI);

    return wt;
}

double ThreeBodyMapper::Isotropic2Weight(const FourVector &p1, const FourVector &p2, double &ran1,
                                         double &ran2, double ctmin, double ctmax) {
    FourVector p1h, p = p1 + p2;

    Boost(1, p, p1, p1h);
    ran1 = (p1h.Z() / p1h.P() - ctmin) / (ctmax - ctmin);
    ran2 = asin(p1h.X() / p1h.Pt()) / (2. * M_PI);
    if(p1h.Y() < 0.) ran2 = .5 - ran2;
    if(ran2 < 0.) ran2 += 1.;

    double massfactor = SqLam(p.M2(), p1.M2(), p2.M2());
    if(massfactor < 1.e-12) return 0.;
    return 2. / M_PI / massfactor * 2.0 / (ctmax - ctmin);
}

void ThreeBodyMapper::Boost(int lflag, const FourVector &q, const FourVector &ph, FourVector &p) {
    if(q.M2() < 0.) { return; }
    double rsq = sqrt(q.M2());
    if(lflag == 0) {
        p.SetE((q.E() * ph.E() + q.Vec3() * ph.Vec3()) / rsq);
        double c1 = (ph.E() + p.E()) / (rsq + q.E());
        p = FourVector(ph.Vec3() + c1 * q.Vec3(), p.E());
    } else {
        p.SetE(q * ph / rsq);
        double c1 = (p.E() + ph.E()) / (rsq + q.E());
        p = FourVector(ph.Vec3() - c1 * q.Vec3(), p.E());
    }
}

#ifdef ACHILLES_SHERPA_INTERFACE
void SherpaMapper::GeneratePoint(std::vector<FourVector> &point, const std::vector<double> &rans) {
    std::vector<Vec4D> mom(point.size());
    mom[0] =
        Vec4D(point[0][0] / 1_GeV, point[0][1] / 1_GeV, point[0][2] / 1_GeV, point[0][3] / 1_GeV);
    mom[1] =
        Vec4D(point[1][0] / 1_GeV, point[1][1] / 1_GeV, point[1][2] / 1_GeV, point[1][3] / 1_GeV);
    sherpa_mapper->GeneratePoint(mom, rans);
    for(size_t i = 2; i < point.size(); ++i) {
        point[i] =
            FourVector(mom[i][0] * 1_GeV, mom[i][1] * 1_GeV, mom[i][2] * 1_GeV, mom[i][3] * 1_GeV);
    }
}

double SherpaMapper::GenerateWeight(const std::vector<FourVector> &point,
                                    std::vector<double> &rans) {
    std::vector<Vec4D> mom{};
    for(const auto &pt : point)
        mom.emplace_back(pt[0] / 1_GeV, pt[1] / 1_GeV, pt[2] / 1_GeV, pt[3] / 1_GeV);
    return sherpa_mapper->GenerateWeight(mom, rans);
}
#endif
