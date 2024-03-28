#ifndef MAPPER_HH
#define MAPPER_HH

#include <memory>
#include <string>
#include <vector>

#include "spdlog/spdlog.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "yaml-cpp/yaml.h"
#pragma GCC diagnostic pop

#include "Achilles/YAML/Beams.hh"

namespace achilles {

template <typename T> class Mapper {
  public:
    template <typename C> using Mapper_ptr = std::unique_ptr<Mapper<C>>;
    template <typename C> using Mapper_sptr = std::shared_ptr<Mapper<C>>;

    Mapper() = default;
    Mapper(std::string name) : mapping_name{std::move(name)} {}
    Mapper(const Mapper &) = delete;
    Mapper(Mapper &&) = delete;
    Mapper &operator=(const Mapper &) = delete;
    Mapper &operator=(Mapper &&) = delete;
    virtual ~Mapper() = default;

    // Functions
    virtual void GeneratePoint(std::vector<T> &, const std::vector<double> &) = 0;
    virtual double GenerateWeight(const std::vector<T> &, std::vector<double> &) = 0;
    virtual size_t NDims() const = 0;
    virtual void SetMasses(std::vector<double> masses) { m_masses = std::move(masses); }
    virtual const std::vector<double> &Masses() const { return m_masses; }
    virtual void SetGaugeBosonMass(double) {}

    // Printers
    static void Print(const char *func, const std::vector<T> &point,
                      const std::vector<double> &rans) {
        spdlog::trace("{}", func);
        size_t idx = 0;
        spdlog::trace("  Point:");
        for(const auto &p : point) spdlog::trace("    - {}: {}", idx++, p);
        idx = 0;
        spdlog::trace("  Rans:");
        for(const auto &r : rans) spdlog::trace("    - {}: {}", idx++, r);
    }
    virtual YAML::Node ToYAML() const = 0;

  protected:
    double Smin() const {
        double smin = 0;
        for(const auto &mass : Masses()) { smin += sqrt(mass); }
        return smin * smin;
    }

  private:
    std::string mapping_name{};
    std::vector<double> m_masses;
};

} // namespace achilles

#endif
