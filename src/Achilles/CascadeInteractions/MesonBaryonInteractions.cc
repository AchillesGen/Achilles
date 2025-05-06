#include "Achilles/CascadeInteractions/MesonBaryonInteractions.hh"
#include "Achilles/Event.hh"
#include "Achilles/Poincare.hh"
#include "Achilles/Random.hh"
#include <limits>

using namespace achilles;

InteractionResults MesonBaryonInteraction::CrossSection(Event &event, size_t part1,
                                                        size_t part2) const {
    const auto &particle1 = event.Hadrons()[part1];
    const auto &particle2 = event.Hadrons()[part2];

    int pidm = particle1.ID();
    int pidb = particle2.ID(); // Usually the case

    if(particle1.Info().IntSpin() % 2 == 1) {
        // In this case what we thought was a meson has half-integer spin
        pidb = particle2.ID();
        pidm = particle1.ID();
    }

    // Could encapsulate all the following in one function, but here we make use of achilles
    // utilities
    size_t ichan = m_amps.GetCchannel(pidm, pidb);

    InteractionResults results;
    // return empty interaction list, no initial channel matches
    if(ichan == std::numeric_limits<size_t>::max()) { return results; }

    double W = (particle1.Momentum() + particle2.Momentum()).M();

    auto CS_fchan = m_amps.GetAllCSW(ichan, W);

    // Remove cross sections that are identically zero
    for(auto &CS_i : CS_fchan) {
        double CS = CS_i.first;
        if(CS > 0.) {
            size_t fchan = CS_i.second;
            int meson_out = m_amps.MesonPID_Cchan(fchan);
            int baryon_out = m_amps.BaryonPID_Cchan(fchan);
            results.push_back({{PID(meson_out), PID(baryon_out)}, CS});
        }
    }

    return results;
}

std::vector<std::pair<PID, PID>> MesonBaryonInteraction::InitialStates() const {
    std::vector<std::pair<PID, PID>> states;

    size_t nChan = m_amps.NChargeChannels();

    for(size_t i = 0; i < nChan; i++) {
        int pidm = m_amps.MesonPID_Cchan(i);
        int pidb = m_amps.BaryonPID_Cchan(i);
        states.push_back({PID(pidm), PID(pidb)});
    }
    return states;
}

std::vector<Particle> MesonBaryonInteraction::GenerateMomentum(const Particle &particle1,
                                                               const Particle &particle2,
                                                               const std::vector<PID> &out_pids,
                                                               Random &random) const {
    // PIDs from initial states
    int ipidm = particle1.ID();
    int ipidb = particle2.ID(); // Usually the case

    if(particle1.Info().IntSpin() % 2 == 1) {
        // In this case what we thought was a meson has half-integer spin
        ipidb = particle2.ID();
        ipidm = particle1.ID();
    }

    size_t ichan = m_amps.GetCchannel(ipidm, ipidb); // channel in

    int fpidm = int(out_pids[0]);
    int fpidb = int(out_pids[1]);

    if(fpidm % 2 == 0) // Works for all, last number in pid should be 2*j except for the KS and KL ?
                       // these are not included here anyway
    {
        // In this case what we thought was a meson has half-integer spin
        fpidb = fpidm;
        fpidm = int(out_pids[1]);
    }

    // Channel out
    size_t fchan = m_amps.GetCchannel(fpidm, fpidb);

    double W = (particle1.Momentum() + particle2.Momentum()).M();

    // Get angular distribution
    f_Polynomial poly_angles = m_amps.Get_CSpoly_W(W, ichan, fchan);

    // We convert the CS to the cdf later, could in principle get the CDF directly
    // For this we need the values:
    double CSmin = poly_angles.eval_If(-1.);
    double CStot = poly_angles.eval_If(1.) - CSmin;

    std::vector<double> rans(3);
    random.Generate(rans); // Between 0 and 1 ?

    double x = rans[0];
    //    if (x < 0 || x > 1) BIG Problem

    // Bisection method, should be incapsulated elsewhere in the future
    //
    // Parameters for bisection
    int Nmax = 40; // too much : Nmax*tol > 1 => should converge
    int Nit = 0;
    double tol = 0.025;   // tolerance on cosine
    double eps = 0.00001; // tolerance on f(x) - c = 0

    double fc, c;
    double b = 1;
    // random point to start bisection and get random bin definitions:
    double a = -1 + rans[2] * 2;
    double fa = (poly_angles.eval_If(a) - CSmin) / CStot - x;

    if(std::abs(fa) < tol) {
        // Found the solution
        c = a;
        Nit = Nmax + 1;
    } else if(fa > 0.) {
        // a is the upper limit instead of the lower
        b = a;
        a = -1;
        fa = -x;
    }

    while(Nit < Nmax) {
        c = (a + b) / 2.;
        fc = (poly_angles.eval_If(c) - CSmin) / CStot - x;
        if((std::abs(fc) < eps) || ((b - a) / 2. < tol)) { break; }
        if(fc / fa > 0) {
            a = c;
            fa = fc;
        } else {
            b = c;
        }
        Nit++;
    }

    // END bisection

    double cos_CMS = c;
    double phi_CMS = rans[1] * 2 * M_PI;
    double sin_CMS = sqrt(1. - c * c);
    double sinphi = sin(phi_CMS);
    double cosphi = cos(phi_CMS);

    // Could generate state here directly if outgoing state is same as initial
    // More generally below, for arbitrary meson-baryon system created

    // Masses of final-state particles
    // Use particleInfo
    double mM = ParticleInfo(fpidm).Mass();
    double mB = ParticleInfo(fpidb).Mass();

    // CMS momenta and energy of final-state particles
    // double s = W * W;
    // double mM2 = mM * mM;
    // double mB2 = mB * mB;
    // double pfCMS2 =
    //     1. / 4 / s * (s * s + mM2 * mM2 + mB2 * mB2 - 2 * s * (mM2 + mB2) - 2 * mM2 * mB2);
    // double pfCMS = sqrt(pfCMS2);
    auto p01 = (particle1.Momentum() + particle2.Momentum());
    auto s = p01.M2();
    auto sqrts = sqrt(s);
    auto boostVec = p01.BoostVector();
    auto mom0 = particle1.Momentum().Boost(-boostVec);
    Poincare zax(mom0, FourVector(1., 0., 0., 1.));

    double EmCMS = sqrts / 2 * (1 + mM * mM / s - mB * mB / s);
    double EbCMS = sqrts / 2 * (1 - mM * mM / s + mB * mB / s);
    auto lambda = sqrt(pow(s - mM * mM - mB * mB, 2) - 4 * mM * mM * mB * mB);
    auto pfCMS = lambda / 2 / sqrts;

    // Generate outgoing momentum
    FourVector p1Out =
        FourVector(EmCMS, pfCMS * sin_CMS * cosphi, pfCMS * sin_CMS * sinphi, pfCMS * cos_CMS);
    FourVector p2Out =
        FourVector(EbCMS, -pfCMS * sin_CMS * cosphi, -pfCMS * sin_CMS * sinphi, -pfCMS * cos_CMS);

    // Rotate back
    zax.RotateBack(p1Out);
    zax.RotateBack(p2Out);

    // // Boost
    // ThreeVector boostCM = (particle1.Momentum() + particle2.Momentum()).BoostVector();
    // FourVector p1Lab = particle1.Momentum();
    // FourVector p1CM = p1Lab.Boost(-boostCM);

    // // The angle is defined with respect to pCMS_initial
    // // We construct an orthogonal coordinate system
    // ThreeVector p_axis = p1CM.Vec3().Unit();
    // ThreeVector p_perp;

    // if(p_axis.Py() == 0.) {
    //     // An arbitrary orthogonal vector is y
    //     p_perp = {0., 1., 0.};
    // } else if(p_axis.Px() == 0) {
    //     // An arbitrary orthogonal vector is x
    //     p_perp = {1., 0., 0.};
    // } else {
    //     double norm = sqrt(1. - p_axis.Pz() * p_axis.Pz());
    //     // An arbitrary orthonormal vector is
    //     p_perp = {p_axis[1] / norm, -p_axis[0] / norm, 0.};
    // }

    // ThreeVector p_out = pfCMS * (p_axis * cos_CMS + p_perp * sin_CMS * sinphi +
    //                              p_perp.Cross(p_axis) * sin_CMS * cosphi); // Rotate p_out

    // Boost back to lab frame
    p1Out = p1Out.Boost(boostVec);
    p2Out = p2Out.Boost(boostVec);

    // Not sure what to do with the positions here? !!! ?
    return {Particle{out_pids[0], p1Out, particle1.Position()},
            Particle{out_pids[1], p2Out, particle2.Position()}};
}
