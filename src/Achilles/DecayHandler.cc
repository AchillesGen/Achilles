#include "Achilles/DecayHandler.hh"
#include "Achilles/MesonBaryonAmplitudes.hh"
#include "Achilles/Particle.hh"
#include "Achilles/Poincare.hh"
#include "Achilles/Random.hh"

#include "spdlog/spdlog.h"

using achilles::DecayHandler;

std::ofstream DecayHandler::m_angular_dist_file{"angular_distribution.txt"};

DecayHandler::DecayHandler(const std::string &filename, double tolerance) {
    auto node = YAML::LoadFile(filename);

    for(const auto &particle : node) {
        auto pid = particle.first.as<PID>();
        double br_total = 0;
        for(const auto &decay_mode : particle.second) {
            DecayMode decay{decay_mode["BranchingRatio"].as<double>(),
                            decay_mode["AngularMom"].as<size_t>(),
                            decay_mode["OutIDs"].as<std::vector<PID>>()};
            br_total += decay.branching_ratio;
            std::sort(decay.out_ids.begin(), decay.out_ids.end());
            m_decays[pid].push_back(decay);
        }

        if(std::abs(br_total - 1) > tolerance) {
            spdlog::warn(
                "DecayHandler: Branching ratio not within {} of 1 ({}) for PID {}. Rescaling.",
                tolerance, br_total, pid);
            for(auto &mode : m_decays[pid]) { mode.branching_ratio /= br_total; }
        }
    }
}

std::vector<achilles::Particle> DecayHandler::Decay(const Particle &part) const {
    spdlog::debug("Decaying: {}", part);
    auto boost = part.Momentum().BoostVector();
    double m2 = part.Momentum().M2();

    // Select decay channel
    auto ratios = BranchingRatios(part);
    spdlog::trace("Branching Rations: {}", fmt::join(ratios, ","));
    // Remove kinematically disallowed decays
    size_t tmp = 0;
    for(const auto &mode : m_decays.at(part.ID())) {
        if(mode.out_ids.size() == 2) {
            spdlog::trace("Checking decay {} ({}) -> {} ({}) + {} ({})", part.ID(), sqrt(m2),
                          mode.out_ids[0], ParticleInfo(mode.out_ids[0]).Mass(), mode.out_ids[1],
                          ParticleInfo(mode.out_ids[1]).Mass());
            auto mass1 = ParticleInfo(mode.out_ids[0]).Mass();
            auto mass2 = ParticleInfo(mode.out_ids[1]).Mass();
            if(mass1 + mass2 > sqrt(m2)) {
                spdlog::debug("Removing decay {} ({}) -> {} ({}) + {} ({})", part.ID(), sqrt(m2),
                              mode.out_ids[0], mass1, mode.out_ids[1], mass2);
                ratios[tmp] = 0;
            }
        }
        tmp++;
    }
    auto idx = Random::Instance().SelectIndex(ratios);
    auto mode = m_decays.at(part.ID())[idx];
    spdlog::trace("Decaying to mode: {}", mode);

    std::vector<Particle> outgoing;
    if(mode.out_ids.size() == 2) {
        if(part.Mothers().size() != 2 || mode.out_ids[0] == PID::photon() ||
           mode.out_ids[1] == PID::photon()) {
            outgoing = TwoBodyDecay(m2, mode.out_ids, mode.angular_mom);
        } else {
            static MBAmplitudes m_amps;
            auto particle1 = part.Mothers()[0];
            auto particle2 = part.Mothers()[1];
            // PIDs from initial states
            int ipidm = particle1.ID();
            int ipidb = particle2.ID(); // Usually the case

            if(particle1.Info().IntSpin() % 2 == 1) {
                // In this case what we thought was a meson has half-integer spin
                ipidb = particle2.ID();
                ipidm = particle1.ID();
            }

            size_t ichan = m_amps.GetCchannel(ipidm, ipidb); // channel in

            int fpidm = int(mode.out_ids[0]);
            int fpidb = int(mode.out_ids[1]);

            if(fpidm % 2 == 0) // Works for all, last number in pid should be 2*j except for the KS
                               // and KL ? these are not included here anyway
            {
                // In this case what we thought was a meson has half-integer spin
                fpidb = fpidm;
                fpidm = int(mode.out_ids[1]);
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
            Random::Instance().Generate(rans);

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

            m_angular_dist_file << c << "\n";

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
            //     1. / 4 / s * (s * s + mM2 * mM2 + mB2 * mB2 - 2 * s * (mM2 + mB2) - 2 * mM2 *
            //     mB2);
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
            FourVector p1Out = FourVector(EmCMS, pfCMS * sin_CMS * cosphi, pfCMS * sin_CMS * sinphi,
                                          pfCMS * cos_CMS);
            FourVector p2Out = FourVector(EbCMS, -pfCMS * sin_CMS * cosphi,
                                          -pfCMS * sin_CMS * sinphi, -pfCMS * cos_CMS);

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
            //                              p_perp.Cross(p_axis) * sin_CMS * cosphi); // Rotate
            //                              p_out

            // Boost back to lab frame
            p1Out = p1Out.Boost(boostVec);
            p2Out = p2Out.Boost(boostVec);

            // Not sure what to do with the positions here? !!! ?
            return {Particle{mode.out_ids[0], p1Out, part.Position()},
                    Particle{mode.out_ids[1], p2Out, part.Position()}};
        }
    } else {
        throw std::runtime_error("Three body and above decays are not implemented yet");
    }

    // Rotate to have z-axis along the decaying particle's momentum direction
    if(!part.Mothers().empty()) {
        spdlog::debug("Mother[0]: {}, Momentum: {}", part.Mothers()[0].ID(),
                      part.Mothers()[0].Momentum());
        auto mom0 = part.Mothers()[0].Momentum().Boost(-boost);
        spdlog::debug("Boosted Mom[0]: {}", mom0);
        Poincare zax(mom0, FourVector{1.0, 0.0, 0.0, 1.0});
        for(auto &particle : outgoing) {
            spdlog::debug("Before Rotation: {}, Momentum: {}", particle.ID(), particle.Momentum());
            zax.RotateBack(particle.Momentum());
            spdlog::debug("After Rotation: {}, Momentum: {} ({})", particle.ID(),
                          particle.Momentum(), particle.Momentum().P());
        }
    }

    static double total_mom = 0;
    static size_t count = 0;
    // Boost back to lab frame
    for(auto &particle : outgoing) {
        particle.Momentum() = particle.Momentum().Boost(boost);
        spdlog::debug("After Boost: {}, {} ({})", particle.ID(), particle.Momentum(),
                      particle.Momentum().P());
        if(particle.Info().IsPion()) {
            total_mom += particle.Momentum().P();
            count++;
            if(count % 1000 == 0) {
                spdlog::info("Average pion momentum: {}", total_mom / static_cast<double>(count));
            }
        }
    }
    return outgoing;
}

std::vector<double> DecayHandler::BranchingRatios(const Particle &part) const {
    std::vector<double> ratios;
    for(const auto &mode : m_decays.at(part.ID())) { ratios.push_back(mode.branching_ratio); }
    return ratios;
}

std::vector<achilles::DecayMode> DecayHandler::AllowedDecays(PID pid) const {
    return m_decays.at(pid);
}

double DecayHandler::BranchingRatio(PID res, std::vector<PID> out) const {
    spdlog::trace("Getting BranchingRatio for {} -> {}", res, fmt::join(out, ","));
    const auto &decays = m_decays.at(res);
    std::sort(out.begin(), out.end());
    for(const auto &decay : decays) {
        if(decay.out_ids == out) return decay.branching_ratio;
    }
    return 0;
}

std::vector<achilles::Particle>
DecayHandler::TwoBodyDecay(double mass2, const std::vector<PID> &pids, size_t angular_mom) const {
    // TODO: Add in angular momentum (i.e. learn more about how Sherpa handles this)
    double sqrts = sqrt(mass2);

    auto ma = ParticleInfo(pids[0]).Mass();
    auto mb = ParticleInfo(pids[1]).Mass();

    const double Ea = sqrts / 2. * (1 + ma * ma / mass2 - mb * mb / mass2);
    const double Eb = sqrts / 2. * (1 + mb * mb / mass2 - ma * ma / mass2);
    auto lambda = sqrt(pow(mass2 - ma * ma - mb * mb, 2) - 4 * ma * ma * mb * mb);
    auto pf = lambda / 2 / sqrts;

    std::vector<double> rans(2);
    Random::Instance().Generate(rans);

    double cost;

    if(angular_mom == 2) {
        double term = cbrt(9 - 18 * rans[0] +
                           2 * sqrt(3.0) * sqrt(7.0 - 27 * rans[0] + 27 * rans[0] * rans[0]));
        cost = 1.0 / (cbrt(3.0) * term) - term / cbrt(9);
    } else {
        cost = 2 * rans[0] - 1;
    }

    double sint = sqrt(1 - cost * cost);
    double phi = 2 * M_PI * rans[1];

    Particle part1{pids[0], {Ea, pf * sint * cos(phi), pf * sint * sin(phi), pf * cost}};
    Particle part2{pids[1], {Eb, -pf * sint * cos(phi), -pf * sint * sin(phi), -pf * cost}};

    return {part1, part2};
}
