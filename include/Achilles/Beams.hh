#ifndef BEAMS_HH
#define BEAMS_HH

#include "Achilles/Achilles.hh"
#include "Achilles/Histogram.hh"
#include <memory>
#include <set>

#include "Achilles/FourVector.hh"
#include "Achilles/PDFBase.hh"
#include "Achilles/ParticleInfo.hh"

namespace achilles {

class FluxType {
  public:
    FluxType() = default;
    FluxType(const FluxType &) = default;
    FluxType(FluxType &&) = default;
    FluxType &operator=(const FluxType &) = default;
    FluxType &operator=(FluxType &&) = default;
    virtual ~FluxType() = default;

    virtual int NVariables() const = 0;
    virtual FourVector Flux(const std::vector<double> &, double) const = 0;
    virtual double GenerateWeight(const FourVector &, std::vector<double> &, double) const = 0;
    virtual std::string Type() const = 0;
    virtual double MinEnergy() const = 0;
    virtual double MaxEnergy() const = 0;
    virtual double EvaluateFlux(const FourVector &) const = 0;
};

class Monochromatic : public FluxType {
  public:
    Monochromatic(const double &energy) : m_energy(energy) {}
    int NVariables() const override { return 0; }
    FourVector Flux(const std::vector<double> &, double) const override {
        return {m_energy, 0, 0, m_energy};
    }
    double GenerateWeight(const FourVector &, std::vector<double> &, double) const override {
        return 1;
    }
    std::string Type() const override { return "Monochromatic"; }
    double MinEnergy() const override { return m_energy; }
    double MaxEnergy() const override { return m_energy; }
    double EvaluateFlux(const FourVector &) const override { return 1; }

  private:
    double m_energy;
};

// TODO: Figure out how to handle spectrums
// How to handle this correctly?
// 1. Generate the energy according to some distribution
// 2. A. Multiply maximum energy by a fraction to get neutrino energy
//    B. Calculate the probability for this energy
// The first requires being able to appropriately sample from the distribution
// (i.e. having a form for the inverse of the CDF). The second would take advantage
// of the importance sampling of Vegas. Furthermore, the second would make it easier to
// combine beams of different initial state particles in a straightforward manner
class Spectrum : public FluxType {
  public:
    enum class Type {
        Histogram,
        ROOTHist,
    };
    enum class FluxFormat {
        Achilles,
        MiniBooNE,
        T2K,
    };

    Spectrum(const YAML::Node &);
    int NVariables() const override { return 1; }
    FourVector Flux(const std::vector<double> &, double) const override;
    double GenerateWeight(const FourVector &, std::vector<double> &, double) const override;
    std::string Type() const override { return "Spectrum"; }
    std::string Format() const;
    double MinEnergy() const override { return m_min_energy; }
    double MaxEnergy() const override { return m_max_energy; }
    double EvaluateFlux(const FourVector &) const override;

  private:
    void AchillesHeader(std::ifstream &);
    void MiniBooNEHeader(std::ifstream &);
    void T2KHeader(std::ifstream &);

    enum class flux_units {
        v_m2_POT_500MeV,
        v_nb_POT_MeV,
        v_cm2_POT_MeV,
        v_cm2_POT_50MeV,
        cm2_50MeV,
    };
    std::function<double(double)> m_flux{};
    double m_min_energy{}, m_max_energy{};
    double m_delta_energy{}, m_energy_units{1};
    double m_flux_integral{};
    flux_units m_units;
    FluxFormat m_format;
};

class PDFBeam : public FluxType {
  public:
    PDFBeam(const YAML::Node &);
    int NVariables() const override { return 1; }
    FourVector Flux(const std::vector<double> &, double) const override;
    double GenerateWeight(const FourVector &, std::vector<double> &, double) const override;
    std::string Type() const override { return "PDFBeam"; }
    double MinEnergy() const override { return 0; }
    double MaxEnergy() const override { return 0; }
    double EvaluateFlux(const FourVector &) const override;

  private:
    std::unique_ptr<PDFBase> p_pdf;
};

class Beam {
  public:
    using BeamMap = std::map<achilles::PID, std::shared_ptr<FluxType>>;

    Beam(BeamMap beams);
    Beam(const Beam &) = delete;
    Beam(Beam &&) = default;
    Beam &operator=(const Beam &) = delete;
    Beam &operator=(Beam &&) = default;
    virtual ~Beam() = default;

    Beam() { n_vars = 0; }
    MOCK int NVariables() const { return n_vars; }
    MOCK FourVector Flux(const PID pid, const std::vector<double> &rans, double smin) const {
        return m_beams.at(pid)->Flux(rans, smin);
    }
    MOCK double GenerateWeight(const PID pid, const FourVector &p, std::vector<double> &rans,
                               double smin) const {
        return m_beams.at(pid)->GenerateWeight(p, rans, smin);
    }
    size_t NBeams() const { return m_beams.size(); }
    MOCK const std::set<PID> &BeamIDs() const { return m_pids; }
    double MaxEnergy() const;
    MOCK double EvaluateFlux(const PID pid, const FourVector &p) const {
        return m_beams.at(pid)->EvaluateFlux(p);
    }

    // Accessors
    std::shared_ptr<FluxType> operator[](const PID pid) { return m_beams[pid]; }
    std::shared_ptr<FluxType> at(const PID pid) const { return m_beams.at(pid); }
    std::shared_ptr<FluxType> operator[](const PID pid) const { return m_beams.at(pid); }

    friend YAML::convert<Beam>;

  private:
    int n_vars;
    std::set<PID> m_pids;
    BeamMap m_beams;
};

} // namespace achilles

#endif
