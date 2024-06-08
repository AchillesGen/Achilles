#include <fstream>
#include <iostream>
#include <map>

#include "Achilles/Potential.hh"
#include "spdlog/spdlog.h"

#include "Achilles/Constants.hh"
#include "Achilles/FourVector.hh"
#include "Achilles/Interactions.hh"
#include "Achilles/Particle.hh"
#include "Achilles/Random.hh"
#include "Achilles/ThreeVector.hh"
#include "Achilles/Utilities.hh"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#ifndef __clang__
#pragma GCC diagnostic ignored "-Wuseless-cast"
#endif
#pragma GCC diagnostic ignored "-Wnull-dereference"
#include "highfive/H5File.hpp"
#pragma GCC diagnostic pop

using namespace achilles;

const std::map<std::string, double> HZETRN = {{"a", 5.0_MeV},
                                              {"b", 0.199 / sqrt(1_MeV)},
                                              {"c", 0.451 * pow(1_MeV, -0.258)},
                                              {"d", 25.0_MeV},
                                              {"e", 134.0_MeV},
                                              {"f", 1.187 * pow(1_MeV, -0.35)},
                                              {"g", 0.1_MeV},
                                              {"h", 0.282_MeV}};

const std::map<std::string, double> PDG = {
    {"Zpp", 33.45_mb},        {"Zpn", 35.80_mb},  {"Y1pp", 42.53_mb}, {"Y1pn", 40.15_mb},
    {"Y2pp", 33.34_mb},       {"Y2pn", 30.00_mb}, {"B", 0.308_mb},    {"s1", 1.0 * pow(1_GeV, 2)},
    {"s0", pow(5.38_GeV, 2)}, {"n1", 0.458},      {"n2", 0.545}};

const std::map<std::string, double> JWN = {{"gamma", 52.5 * pow(1_GeV, 0.16)}, // mb
                                           {"alpha", 0.00369 / 1_MeV},
                                           {"beta", 0.00895741 * pow(1_MeV, -0.8)}};

double Interaction::TotalCrossSection(const Particle &particle1, const Particle &particle2) const {
    auto cross_sections = CrossSection(particle1, particle2);
    double total = 0;
    for(const auto &cross_section : cross_sections) total += cross_section.cross_section;
    return total;
}

double Interaction::CrossSectionLab(bool samePID, const double &pLab) const noexcept {
    const double tLab = sqrt(pow(pLab, 2) + pow(Constant::mN, 2)) - Constant::mN;
    if(samePID) {
        if(pLab < 1.8_GeV) {
            if(tLab >= 25_MeV)
                return (1.0 + HZETRN.at("a") / tLab) *
                       (40 + 109.0 * std::cos(HZETRN.at("b") * sqrt(tLab)) *
                                 exp(-HZETRN.at("c") * pow(tLab - HZETRN.at("d"), 0.258)));
            else
                return exp(6.51 * exp(-pow(tLab / HZETRN.at("e"), 0.7)));
        } else if(pLab <= 4.7_GeV) {
            return JWN.at("gamma") / pow(pLab, 0.16);
        } else {
            double ecm2 =
                2 * Constant::mN * (Constant::mN + sqrt(pow(pLab, 2) + pow(Constant::mN, 2)));
            return PDG.at("Zpp") + PDG.at("B") * pow(log(ecm2 / PDG.at("s0")), 2) +
                   PDG.at("Y1pp") * pow(PDG.at("s1") / ecm2, PDG.at("n1")) -
                   PDG.at("Y2pp") * pow(PDG.at("s1") / ecm2, PDG.at("n2"));
        }
    } else {
        if(pLab < 0.5_GeV) {
            if(tLab >= 0.1_MeV)
                return 38.0 + 12500.0 * exp(-HZETRN.at("f") * pow(tLab - HZETRN.at("g"), 0.35));
            else
                return 26000 * exp(-pow(tLab / HZETRN.at("h"), 0.3));
        } else if(pLab <= 2.0_GeV) {
            return 40 + 10 * cos(JWN.at("alpha") * pLab - 0.943) *
                            exp(-JWN.at("beta") * pow(pLab, 0.8) + 2);
        } else {
            double ecm2 =
                2 * Constant::mN * (Constant::mN + sqrt(pow(pLab, 2) + pow(Constant::mN, 2)));
            return PDG.at("Zpn") + PDG.at("B") * pow(log(ecm2 / PDG.at("s0")), 2) +
                   PDG.at("Y1pn") * pow(PDG.at("s1") / ecm2, PDG.at("n1")) -
                   PDG.at("Y2pn") * pow(PDG.at("s1") / ecm2, PDG.at("n2"));
        }
    }
}

GeantInteraction::GeantInteraction(const YAML::Node &node) {
    auto filename = node["GeantData"].as<std::string>();

    // Initialize theta vector
    constexpr double thetaMin = 0.5;
    constexpr double thetaMax = 179.5;
    constexpr size_t nTheta = 180;
    m_theta = Linspace(thetaMin, thetaMax, nTheta);

    constexpr double cdfMin = -3;
    constexpr double cdfMax = 0;
    constexpr size_t nCDF = 200;
    m_cdf = Logspace(cdfMin, cdfMax, nCDF);

    // Read in the Geant4 hdf5 file and get the np and pp groups
    spdlog::info("GeantInteractions: Loading Geant4 data from {0}.", filename);
    HighFive::File file(filename, HighFive::File::ReadOnly);
    HighFive::Group dataNP(file.getGroup("np"));
    HighFive::Group dataPP(file.getGroup("pp"));

    // Get the datasets for np and load into local variables
    LoadData(false, dataNP);

    // Get the datasets for pp and load into local variables
    LoadData(true, dataPP);
}

void GeantInteraction::LoadData(bool samePID, const HighFive::Group &group) {
    // Load datasets
    HighFive::DataSet pcm(group.getDataSet("pcm"));
    HighFive::DataSet sigTot(group.getDataSet("sigtot"));
    HighFive::DataSet sig(group.getDataSet("sig"));

    // Get data for center of momentum
    std::vector<double> pcmVec;
    pcm.read(pcmVec);

    // Get data for total cross-section
    std::vector<double> sigTotVec;
    sigTot.read(sigTotVec);

    // Get data for angular cross-section
    auto dims = sig.getDimensions();
    std::vector<double> sigAngular(dims[0] * dims[1]);
    sig.read(sigAngular.data());

    // Perform interpolation for angles
    achilles::Interp2D interp(pcmVec, m_theta, sigAngular);
    interp.BicubicSpline();

    std::vector<double> theta(pcmVec.size() * m_cdf.size());
    constexpr double accuracy = 1E-6;
    for(size_t i = 0; i < pcmVec.size(); ++i) {
        for(size_t j = 0; j < m_cdf.size(); ++j) {
            auto func = [interp, pcmVec, this, i, j](double x) {
                return interp(pcmVec[i], x) - this->m_cdf[j];
            };
            achilles::Brent brent(func, accuracy);
            if(j != m_cdf.size() - 1) try {
                    theta[i * m_cdf.size() + j] = brent.CalcRoot(m_theta.front(), m_theta.back());
                } catch(std::domain_error &e) {
                    theta[i * m_cdf.size() + j] =
                        m_theta.front() / sigAngular[i * 180 + j] * m_cdf[j];
                }
            else
                theta[i * m_cdf.size() + j] = m_theta.back();
        }
    }

    if(samePID) {
        m_pcmPP = pcmVec;
        m_xsecPP = sigTotVec;
        m_crossSectionPP.SetData(pcmVec, sigTotVec);
        m_crossSectionPP.CubicSpline();
        m_thetaDistPP.SetData(pcmVec, m_cdf, theta);
        m_thetaDistPP.BicubicSpline();
    } else {
        m_pcmNP = pcmVec;
        m_xsecNP = sigTotVec;
        m_crossSectionNP.SetData(pcmVec, sigTotVec);
        m_crossSectionNP.CubicSpline();
        m_thetaDistNP.SetData(pcmVec, m_cdf, theta);
        m_thetaDistNP.BicubicSpline();
    }
}

InteractionResults GeantInteraction::CrossSection(const Particle &particle1,
                                                  const Particle &particle2) const {
    bool samePID = particle1.ID() == particle2.ID();
    ThreeVector boostCM = (particle1.Momentum() + particle2.Momentum()).BoostVector();
    FourVector p1Lab = particle1.Momentum(); //, p2Lab = particle2.Momentum();
    FourVector p1CM = p1Lab.Boost(-boostCM);
    // Generate outgoing momentum
    const double pcm = p1CM.Vec3().Magnitude();

    try {
        if(samePID) {
            return {{{particle1.ID(), particle2.ID()}, m_crossSectionPP(pcm / 1_GeV)}};
        } else {
            return {{{particle1.ID(), particle2.ID()}, m_crossSectionNP(pcm / 1_GeV)}};
        }
    } catch(std::domain_error &e) {
        spdlog::trace("Using Nasa Interaction");
        // double s = (p1Lab+p2Lab).M2();
        double s = (particle1.Momentum() + particle2.Momentum()).M2();
        double smin = pow(particle1.Mass(), 2) + pow(particle2.Mass(), 2);
        double plab = sqrt(pow(s, 2) / smin - s);
        return {{{particle1.ID(), particle2.ID()}, Interaction::CrossSectionLab(samePID, plab)}};
    }
}

std::vector<Particle> GeantInteraction::GenerateMomentum(const Particle &particle1,
                                                         const Particle &particle2,
                                                         const std::vector<PID> &out_pids,
                                                         Random &random) const {
    // Boost to center of mass
    ThreeVector boostCM = (particle1.Momentum() + particle2.Momentum()).BoostVector();
    FourVector p1Lab = particle1.Momentum();
    FourVector p1CM = p1Lab.Boost(-boostCM);

    // Generate outgoing momentum
    bool samePID = particle1.ID() == particle2.ID();
    const double pcm = p1CM.Vec3().Magnitude();
    std::vector<double> rans(2);
    random.Generate(rans);
    ThreeVector momentum = MakeMomentum(samePID, pcm, rans);

    FourVector p1Out = FourVector(p1CM.E(), momentum[0], momentum[1], momentum[2]);
    FourVector p2Out = FourVector(p1CM.E(), -momentum[0], -momentum[1], -momentum[2]);

    // Boost back to lab frame
    p1Out = p1Out.Boost(boostCM);
    p2Out = p2Out.Boost(boostCM);

    return {{out_pids[0], p1Out, particle1.Position()}, {out_pids[1], p2Out, particle2.Position()}};
}

ThreeVector GeantInteraction::MakeMomentum(bool samePID, double pcm,
                                           const std::vector<double> &rans) const {
    double pR = pcm;
    double pTheta = CrossSectionAngle(samePID, pcm / 1_GeV, rans[0]);
    double pPhi = 2 * M_PI * rans[1];

    return ThreeVector(ToCartesian({pR, pTheta, pPhi}));
}

double GeantInteraction::CrossSectionAngle(bool samePID, const double &energy,
                                           const double &ran) const {
    try {
        if(samePID)
            return m_thetaDistPP(energy, ran);
        else
            return m_thetaDistNP(energy, ran);
    } catch(std::domain_error &e) {
        spdlog::trace("Using flat angular distribution");
        return acos(2 * ran - 1);
    }
}

/*
GeantInteractionsDt::GeantInteractionsDt(const YAML::Node& node) {
    auto filename_pp = node["GeantDataPP"].as<std::string>();
    auto filename_np = node["GeantDataNP"].as<std::string>();

    // Read in the Geant4 dsigma/dt files
    spdlog::info("GeantInteractionsDt: Loading Geant4 dsigma/dt_pp data from
{0}.", filename_pp); LoadData(true, filename_pp);

    spdlog::info("GeantInteractionsDt: Loading Geant4 dsigma/dt_np data from
{0}.", filename_np); LoadData(false, filename_np); throw;
}

void GeantInteractionsDt::LoadData(bool samePID, const std::string &filename) {
    std::ifstream data(filename);

    // Read in CDFs
    const size_t nlines = samePID ? 40 : 39;
    auto cdf = ReadBlock(data, nlines);

    // Read in tmin and tmax
    auto tmin = ReadBlock(data, 2);
    auto tmax = ReadBlock(data, 1);

    // Read in pcm
    auto ecm2 = ReadBlock(data, 2);
    for(auto & e : ecm2)
        e = 4*(e*e + Constant::mN*Constant::mN);

    // Read in lab energy and max sigma
    auto elab = ReadBlock(data, 2);
    auto max_sig = ReadBlock(data, 2);

    // Read in total cross-section
    auto xsec = ReadBlock(data, 2);

    data.close();

    // Obtain the PDF from the CDF
    const size_t npts = cdf.size()/tmin.size();
    std::vector<double> pdf(ecm2.size()*npts);
    for(size_t i = 0; i < ecm2.size(); ++i) {
        double tmp = 0;
        for(size_t j = 0; j < npts; ++j) {
            if(j == 0) pdf[i*npts] = cdf[i*npts];
            else pdf[i*npts + j] = cdf[i*npts+j] - cdf[i*npts+j-1];

            tmp += pdf[i*npts + j];
        }
    }

    // Store data in class
    auto trange = Linspace(0, 1, nlines);
    std::reverse(trange.begin(), trange.end());
    if(samePID) {
        m_crossSectionPP.SetData(ecm2, xsec);
        m_crossSectionPP.CubicSpline();
        m_pdfPP.SetData(ecm2, trange, pdf);
        m_pdfPP.BicubicSpline();
    } else {
        m_crossSectionNP.SetData(ecm2, xsec);
        m_crossSectionNP.CubicSpline();
        m_pdfNP.SetData(ecm2, trange, pdf);
        m_pdfNP.BicubicSpline();
    }
}

std::vector<double> GeantInteractionsDt::ReadBlock(std::ifstream &file, size_t
nlines) const { std::string line{}; std::vector<std::string> tokens; for(size_t
i = 0; i < nlines; ++i) { std::getline(file, line); tokenize(line, tokens);
    }
    std::vector<double> results;
    for(const auto& token : tokens)
        results.emplace_back(std::stod(token));

    return results;
}

double GeantInteractionsDt::CrossSection(const Particle &particle1,
                                         const Particle &particle2) const {
    bool samePID = particle1.ID() == particle2.ID();
    const double ecm2 = (particle1.Momentum() +
particle2.Momentum()).M2()/1_GeV/1_GeV;

    try{
        if(samePID) {
            return m_crossSectionPP(ecm2);
        } else {
            return m_crossSectionNP(ecm2);
        }
    } catch(std::domain_error &e) {
        spdlog::warn("ecm2 {} not in table, using 0 for xsec", ecm2);
        return 0;
    }
}

// TODO: Finish implementing this function
achilles::Interactions::MomentumPair GeantInteractionsDt::FinalizeMomentum(const
Particle &particle1, const Particle &particle2, std::shared_ptr<Potential> pot)
const {
    // Generate outgoing momentum
    // bool samePID = particle1.ID() == particle2.ID();
    // const double ecm2 = (particle1.Momentum() + particle2.Momentum()).M2();
    // const double tmin = 2*particle1.Mass() + 2*particle2.Mass() - ecm2;
    auto p1 = particle1.Momentum();
    auto p2 = particle2.Momentum();
    // auto q_free = p1 + p2;
    p1.E() = pot->Hamiltonian(p1.P(), particle1.Radius());
    p2.E() = pot->Hamiltonian(p2.P(), particle2.Radius());
    // auto q = p1+p2;

    // Rotate so (p1 + p2) is along the z-axis
    // auto rotation  = q.AlignZ();
    // q = q.Rotate(rotation);

    // Generate random variables
    // std::array<double, 2> rans{};
    // Random::Instance().Generate(rans, 0.0, 1.0);
    // const double phi = 2*M_PI*rans[0];
    // const double t =  samePID ?
    //     (m_pdfPP(ecm2, rans[1])-1)*tmin : (m_pdfNP(ecm2, rans[1])-1)*tmin;

    // ThreeVector momentum = MakeMomentum(samePID, ecm2, rans);

    return {p1, p2};
}

// TODO: Finish implementing this function
//ThreeVector GeantInteractionsDt::MakeMomentum(bool samePID,
//                                              const double &pcm,
//                                              const std::array<double, 2>
&rans) const {
//    double pR = pcm;
//}
*/

InteractionResults NasaInteraction::CrossSection(const Particle &particle1,
                                                 const Particle &particle2) const {
    bool samePID = particle1.ID() == particle2.ID();
    FourVector p1Lab = particle1.Momentum(), p2Lab = particle2.Momentum();
    // Generate outgoing momentum
    double s = (p1Lab + p2Lab).M2();
    double smin = pow(particle1.Mass() + particle2.Mass(), 2);
    double plab = sqrt(pow(s, 2) / smin - s);
    return {{{particle1.ID(), particle2.ID()}, CrossSectionLab(samePID, plab)}};
}

std::vector<Particle> NasaInteraction::GenerateMomentum(const Particle &particle1,
                                                        const Particle &particle2,
                                                        const std::vector<PID> &out_pids,
                                                        Random &random) const {
    // Boost to center of mass
    ThreeVector boostCM = (particle1.Momentum() + particle2.Momentum()).BoostVector();
    FourVector p1Lab = particle1.Momentum();
    FourVector p1CM = p1Lab.Boost(-boostCM);

    // Generate outgoing momentum
    bool samePID = particle1.ID() == particle2.ID();
    const double pcm = p1CM.Vec3().Magnitude();
    std::vector<double> rans(2);
    random.Generate(rans);
    ThreeVector momentum = MakeMomentum(samePID, pcm, rans);

    FourVector p1Out = FourVector(p1CM.E(), momentum[0], momentum[1], momentum[2]);
    FourVector p2Out = FourVector(p1CM.E(), -momentum[0], -momentum[1], -momentum[2]);

    // Boost back to lab frame
    p1Out = p1Out.Boost(boostCM);
    p2Out = p2Out.Boost(boostCM);

    return {{out_pids[0], p1Out, particle1.Position()}, {out_pids[1], p2Out, particle2.Position()}};
}

ThreeVector NasaInteraction::MakeMomentum(bool, double pcm, const std::vector<double> &rans) const {
    double pR = pcm;
    double pTheta = acos(2 * rans[0] - 1);
    double pPhi = 2 * M_PI * rans[1];

    return ThreeVector(ToCartesian({pR, pTheta, pPhi}));
}

ConstantInteraction::ConstantInteraction(const YAML::Node &node) {
    for(const auto &interactions : node["InitialStates"]) {
        auto incoming = std::make_pair<PID, PID>(interactions["Incoming"][0].as<PID>(),
                                                 interactions["Incoming"][1].as<PID>());
        AddInteraction(incoming, interactions["Outgoing"].as<InteractionResults>());
    }
}

InteractionResults ConstantInteraction::CrossSection(const Particle &p1, const Particle &p2) const {
    return m_interactions.at({p1.ID(), p2.ID()});
}

// TODO: Implement generic n-body final state phase space
std::vector<Particle> ConstantInteraction::GenerateMomentum(const Particle &particle1,
                                                            const Particle &particle2,
                                                            const std::vector<PID> &out_pids,
                                                            Random &random) const {
    // Boost to center of mass
    ThreeVector boostCM = (particle1.Momentum() + particle2.Momentum()).BoostVector();
    FourVector p1Lab = particle1.Momentum();
    FourVector p1CM = p1Lab.Boost(-boostCM);

    // Generate outgoing momentum
    bool samePID = particle1.ID() == particle2.ID();
    const double pcm = p1CM.Vec3().Magnitude();
    std::vector<double> rans(2);
    random.Generate(rans);
    ThreeVector momentum = MakeMomentum(samePID, pcm, rans);

    FourVector p1Out = FourVector(p1CM.E(), momentum[0], momentum[1], momentum[2]);
    FourVector p2Out = FourVector(p1CM.E(), -momentum[0], -momentum[1], -momentum[2]);

    // Boost back to lab frame
    p1Out = p1Out.Boost(boostCM);
    p2Out = p2Out.Boost(boostCM);

    return {{out_pids[0], p1Out, particle1.Position()}, {out_pids[1], p2Out, particle2.Position()}};
}

ThreeVector ConstantInteraction::MakeMomentum(bool, double pcm,
                                              const std::vector<double> &rans) const {
    double pR = pcm;
    double pTheta = acos(2 * rans[0] - 1);
    double pPhi = 2 * M_PI * rans[1];

    return ThreeVector(ToCartesian({pR, pTheta, pPhi}));
}


InteractionResults MesonBaryonInteraction::CrossSection(const Particle &particle1,
                                                  const Particle &particle2) const {
    int pidm = particle1.ID();
    int pidb = particle2.ID(); //Usually the case

    if (particle1.Info().IntSpin() % 2 == 1)
    {
	//In this case what we thought was a meson has half-integer spin
	pidb = particle2.ID();
	pidm = particle1.ID();
    }

    //Could encapsulate all the following in one function, but here we make use of achilles utilities
    int ichan = Amplitudes.GetCchannel(pidm, pidb);

    InteractionResults results;
    if (ichan < 0){return results;} //return empty, no initial channel matches

    double W =  (particle1.Momentum() + particle2.Momentum()).M();

    std::vector<std::pair<double, int>> CS_fchan = Amplitudes.GetAllCSW(ichan,W);

    //Remove cross sections that are identically zero
    for (auto & CS_i : CS_fchan)
    {
	double CS = CS_i.first;
	if (CS > 0.){
	    int fchan = CS_i.second;
	    int meson_out = Amplitudes.MesonPID_Cchan(fchan);
	    int baryon_out = Amplitudes.BaryonPID_Cchan(fchan);
	    results.push_back( {{PID(meson_out), PID(baryon_out)}, CS} );
	}
    }

    return results;
}


std::vector<Particle> MesonBaryonInteraction::GenerateMomentum(const Particle &particle1,
                                                         const Particle &particle2,
                                                         const std::vector<PID> &out_pids,
                                                         Random &random) const {

    //PIDs from initial states 
    int ipidm = particle1.ID();
    int ipidb = particle2.ID(); //Usually the case

    if (particle1.Info().IntSpin() % 2 == 1)
    {
	//In this case what we thought was a meson has half-integer spin
	ipidb = particle2.ID();
	ipidm = particle1.ID();
    }

    int ichan = Amplitudes.GetCchannel(ipidm, ipidb); //channel in

    int fpidm = int(out_pids[0]);
    int fpidb = int(out_pids[1]);

    if (fpidm % 2 == 0) //Works for all, last number in pid should be 2*j except for the KS and KL ? these are not included here anyway
    {
	//In this case what we thought was a meson has half-integer spin
	fpidb = fpidm;
	fpidm = int(out_pids[1]);
    }
	
    //Channel out
    int fchan = Amplitudes.GetCchannel(fpidm,fpidb);

    double W =  (particle1.Momentum() + particle2.Momentum()).M();

    //Get angular distribution 
    f_Polynomial poly_angles = Amplitudes.Get_CSpoly_W(W, ichan, fchan);

    //We convert the CS to the cdf later, could in principle get the CDF directly
    //For this we need the values:
    double CSmin = poly_angles.eval_If(-1.);
    double CStot = poly_angles.eval_If(1.) - CSmin;

    std::vector<double> rans(3);
    random.Generate(rans); //Between 0 and 1 ? 
    
    double x = rans[0];
//    if (x < 0 || x > 1) BIG Problem

    //Bisection method, should be incapsulated elsewhere in the future
    //
    //Parameters for bisection
    int Nmax = 40; //too much : Nmax*tol > 1 => should converge
    int Nit = 0;
    double tol = 0.025; //tolerance on cosine 
    double eps = 0.00001; //tolerance on f(x) - c = 0


    double fc, c;
    double b = 1; 
    //random point to start bisection and get random bin definitions:
    double a = -1 + rans[2]*2;
    double fa = (poly_angles.eval_If(a) - CSmin)/CStot - x;

    if (std::abs(fa) < tol){
       //Found the solution
       c = a;
       Nit = Nmax + 1; 
    }else if (fa > 0.){
	//a is the upper limit instead of the lower
	b = a;
	a = -1;
	fa = -x;
    }   

    while (Nit < Nmax)
    {
       c = (a + b)/2.;
       fc = (poly_angles.eval_If(c) - CSmin)/CStot - x;
       if ( ( std::abs(fc)  < eps) || ( (b-a)/2. < tol ) ) {break;}
       if (fc/fa > 0) {
           a = c; 
	   fa = fc;
       }else{
           b = c;
       }
       	   Nit++;
    }

    //END bisection
  
    double cos_CMS = c;
    double phi_CMS = rans[1]* 2 * M_PI;
    double sin_CMS = sqrt(1.-c*c);
    double sinphi = sin(phi_CMS);
    double cosphi = cos(phi_CMS);

    //Could generate state here directly if outgoing state is same as initial
    //More generally below, for arbitrary meson-baryon system created

    //Masses of final-state particles
    //Use particleInfo
    double mM = ParticleInfo(fpidm).Mass();
    double mB = ParticleInfo(fpidb).Mass();

    //CMS momenta and energy of final-state particles
    double s = W*W;
    double mM2 = mM*mM;
    double mB2 = mB*mB;
    double pfCMS2 = 1./4/s * ( s*s + mM2*mM2 + mB2*mB2 - 2*s*(mM2 + mB2) -2*mM2*mB2 );
    double pfCMS = sqrt(pfCMS2);

    double EmCMS = sqrt(pfCMS2 + mM2);
    double EbCMS = W - EmCMS;

    // Boost
    ThreeVector boostCM = (particle1.Momentum() + particle2.Momentum()).BoostVector();
    FourVector p1Lab = particle1.Momentum();
    FourVector p1CM = p1Lab.Boost(-boostCM);

    //The angle is defined with respect to pCMS_initial
    //We construct an orthogonal coordinate system
    ThreeVector p_axis = p1CM.Vec3().Unit();
    ThreeVector p_perp;

    if (p_axis.Py() == 0.){
	//An arbitrary orthogonal vector is y
	p_perp = {0.,1.,0.};
    }else if (p_axis.Px() == 0){
	//An arbitrary orthogonal vector is x
	p_perp = {1.,0.,0.};
    }else{
	double norm = sqrt(1. - p_axis.Pz()*p_axis.Pz());
	//An arbitrary orthonormal vector is
	p_perp = {p_axis[1]/norm, -p_axis[0]/norm,0.};
    }

    ThreeVector p_out = pfCMS*(p_axis*cos_CMS + p_perp*sin_CMS*sinphi + p_perp.Cross(p_axis)*sin_CMS*cosphi); //Rotate p_out 

    FourVector p1Out = FourVector(EmCMS, p_out[0], p_out[1], p_out[2]);
    FourVector p2Out = FourVector(EbCMS, -p_out[0], -p_out[1], -p_out[2]);

    // Boost back to lab frame
    p1Out = p1Out.Boost(boostCM);
    p2Out = p2Out.Boost(boostCM);

    //Not sure what to do with the positions here? !!! ?
    return {{out_pids[0], p1Out, particle1.Position()}, {out_pids[1], p2Out, particle2.Position()}};
}
