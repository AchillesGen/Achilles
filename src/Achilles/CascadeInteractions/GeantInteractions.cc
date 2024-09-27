#include "Achilles/CascadeInteractions/GeantInteractions.hh"
#include "Achilles/Event.hh"
#include "Achilles/Random.hh"
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

std::vector<std::pair<PID, PID>> GeantInteraction::InitialStates() const {
    return {{PID::proton(), PID::proton()},
            {PID::neutron(), PID::proton()},
            {PID::neutron(), PID::neutron()}};
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

InteractionResults GeantInteraction::CrossSection(Event &event, size_t part1, size_t part2) const {
    const auto &particle1 = event.Hadrons()[part1];
    const auto &particle2 = event.Hadrons()[part2];

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

    return {Particle{out_pids[0], p1Out, particle1.Position()},
            Particle{out_pids[1], p2Out, particle2.Position()}};
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
