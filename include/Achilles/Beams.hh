#ifndef BEAMS_HH
#define BEAMS_HH

#include "Achilles/Achilles.hh"
#include "Achilles/Factory.hh"
#include "Achilles/Hashable.hh"
#include "Achilles/Histogram.hh"
#include <memory>
#include <set>

#include "Achilles/FourVector.hh"
#include "Achilles/PDFBase.hh"
#include "Achilles/ParticleInfo.hh"

namespace achilles {

struct FluxData;

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

class Spectrum : public FluxType {
  public:
    explicit Spectrum(const FluxData &);
    int NVariables() const override { return 1; }
    FourVector Flux(const std::vector<double> &, double) const override;
    double GenerateWeight(const FourVector &, std::vector<double> &, double) const override;
    std::string Type() const override { return "Spectrum"; }
    std::string Format() const;
    double MinEnergy() const override { return m_min_energy; }
    double MaxEnergy() const override { return m_max_energy; }
    double EvaluateFlux(const FourVector &) const override;

  private:
    std::function<double(double)> m_flux{};
    double m_min_energy{}, m_max_energy{};
    double m_delta_energy{}, m_energy_units{1};
    double m_flux_integral{};
    std::string m_format_str{"Undefined"};
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

class FlatFlux : public FluxType {
  public:
    FlatFlux(const YAML::Node &);
    int NVariables() const override { return 1; }
    FourVector Flux(const std::vector<double> &, double) const override;
    double GenerateWeight(const FourVector &, std::vector<double> &, double) const override;
    std::string Type() const override { return "FlatFlux"; }
    double MinEnergy() const override { return m_min_energy; }
    double MaxEnergy() const override { return m_max_energy; }
    double EvaluateFlux(const FourVector &) const override { return 1; }

  private:
    double m_min_energy, m_max_energy;
};

class GeometryBeam : public FluxType {
  public:
    enum class Mode {
        WarmUp,
        Generate,
    };

    GeometryBeam(const YAML::Node &);
    int NVariables() const override { return 1; }
    FourVector Flux(const std::vector<double> &, double) const override;
    double GenerateWeight(const FourVector &, std::vector<double> &, double) const override;
    std::string Type() const override { return "GeometryBeam"; }
    double MinEnergy() const override;
    double MaxEnergy() const override;
    double EvaluateFlux(const FourVector &) const override;
    void SetMode(Mode mode) { m_mode = mode; }
    void SetInjected(int nu_pid, const FourVector &);

  private:
    Mode m_mode;
};

// Energy support [Min, Max] shared by every species in a beam.
class EnergyDomain {
  public:
    EnergyDomain() = default;
    EnergyDomain(double min, double max) : m_min(min), m_max(max) {}
    double Min() const { return m_min; }
    double Max() const { return m_max; }
    bool Contains(double energy) const { return energy >= m_min && energy <= m_max; }
    bool operator==(const EnergyDomain &other) const {
        return m_min == other.m_min && m_max == other.m_max;
    }

  private:
    double m_min{};
    double m_max{};
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
    const EnergyDomain &Domain() const { return m_domain; }
    double MinEnergy() const;
    double MaxEnergy() const;
    MOCK double EvaluateFlux(const PID pid, const FourVector &p) const {
        return m_beams.at(pid)->EvaluateFlux(p);
    }

    // Accessors
    std::shared_ptr<FluxType> operator[](const PID pid) { return m_beams[pid]; }
    std::shared_ptr<FluxType> at(const PID pid) const { return m_beams.at(pid); }
    std::shared_ptr<FluxType> operator[](const PID pid) const { return m_beams.at(pid); }

    friend YAML::convert<Beam>;
    friend std::hash<Beam>;

  private:
    int n_vars;
    std::set<PID> m_pids;
    BeamMap m_beams;
    EnergyDomain m_domain;
};

} // namespace achilles

template <> struct std::hash<achilles::FluxType> {
    std::size_t operator()(const achilles::FluxType &b) const {
        size_t seed = 0;
        achilles::utils::hash_combine(seed, b.Type(), b.MinEnergy(), b.MaxEnergy());
        return seed;
    }
};

template <> struct std::hash<achilles::Beam> {
    std::size_t operator()(const achilles::Beam &b) const {
        std::size_t seed = 0;
        for(const auto &beam : b.m_beams) {
            seed ^= std::hash<int>{}(beam.first) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            seed ^= std::hash<achilles::FluxType>{}(*(beam.second.get())) + 0x9e3779b9 +
                    (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

#endif
