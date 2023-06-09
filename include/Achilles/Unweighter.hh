#ifndef UNWEIGHTER_HH
#define UNWEIGHTER_HH

#include "Achilles/Factory.hh"
#include "Achilles/Statistics.hh"

namespace achilles {

class Event;

class Unweighter {
  public:
    Unweighter() = default;
    virtual ~Unweighter() = default;
    Unweighter(const Unweighter &) = default;
    Unweighter(Unweighter &&) = default;
    Unweighter &operator=(const Unweighter &) = default;
    Unweighter &operator=(Unweighter &&) = default;

    virtual void AddEvent(double) = 0;
    virtual double AcceptEvent(double) = 0;
    virtual double MaxValue() = 0;

    double Efficiency() const {
        return static_cast<double>(m_accepted) / static_cast<double>(m_total);
    }
    size_t Accepted() const { return m_accepted; }

  protected:
    size_t m_accepted{}, m_total{};
};

template <typename Derived>
using RegistrableUnweighter = Registrable<Unweighter, Derived, const YAML::Node &>;
using UnweighterFactory = Factory<Unweighter, const YAML::Node &>;

class NoUnweighter : public Unweighter, RegistrableUnweighter<NoUnweighter> {
  public:
    NoUnweighter(const YAML::Node &) {}
    void AddEvent(double weight) override {
        if(weight > m_max_weight) m_max_weight = weight;
    }
    double AcceptEvent(double weight) override {
        m_accepted++;
        m_total++;
        return weight;
    }
    double MaxValue() override { return m_max_weight; }

    // Required factory methods
    static std::unique_ptr<Unweighter> Construct(const YAML::Node &node) {
        return std::make_unique<NoUnweighter>(node);
    }
    static std::string Name() { return "None"; }

  private:
    double m_max_weight{};
};

class PercentileUnweighter : public Unweighter, RegistrableUnweighter<PercentileUnweighter> {
  public:
    PercentileUnweighter(const YAML::Node &);
    void AddEvent(double) override;
    double AcceptEvent(double) override;
    double MaxValue() override { return m_percentile.Get(); }

    // Required factory methods
    static std::unique_ptr<Unweighter> Construct(const YAML::Node &);
    static std::string Name() { return "Percentile"; }

  private:
    Percentile m_percentile;
};

} // namespace achilles

#endif
