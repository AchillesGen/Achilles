// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

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
    virtual void SaveState(std::ostream &) const = 0;
    virtual void LoadState(std::istream &) = 0;

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
    void SaveState(std::ostream &os) const override {
        const auto default_precision{os.precision()};
        os << std::setprecision(std::numeric_limits<double>::max_digits10 + 1);
        os << m_max_weight << " " << m_accepted << " " << m_total << " ";
        os << std::setprecision(static_cast<int>(default_precision));
    }
    void LoadState(std::istream &is) override { is >> m_max_weight >> m_accepted >> m_total; }

    // Required factory methods
    static std::unique_ptr<Unweighter> Construct(const YAML::Node &node) {
        return std::make_unique<NoUnweighter>(node);
    }
    static std::string Name() { return "None"; }

  private:
    double m_max_weight{};
};

class SortedWeightUnweighter : public Unweighter {
  public:
    SortedWeightUnweighter(double param) : m_param{param} {}

    void AddEvent(double weight) override;
    double AcceptEvent(double weight) override;
    double MaxValue() override;
    void SaveState(std::ostream &os) const override;
    void LoadState(std::istream &is) override;

    double RealizedTailFraction();
    double RealizedExcessFraction();
    bool CapAtMaximum();
    bool Frozen() const { return m_frozen; }

  protected:
    virtual double ComputeCap(const std::vector<double> &sorted_weights, double total) const = 0;
    virtual std::string RuleTag() const = 0;
    void EnsureCap();

    std::vector<double> m_weights{};
    double m_param{};
    double m_max_weight{};
    double m_cap{};
    bool m_dirty{true};
    bool m_frozen{false};
};

// Cap = the p-th percentile of |weights|.
class PercentileUnweighter : public SortedWeightUnweighter,
                             RegistrableUnweighter<PercentileUnweighter> {
  public:
    PercentileUnweighter(const YAML::Node &);
    static std::unique_ptr<Unweighter> Construct(const YAML::Node &);
    static std::string Name() { return "Percentile"; }

  protected:
    double ComputeCap(const std::vector<double> &sorted_weights, double total) const override;
    std::string RuleTag() const override { return Name(); }
};

// Cap = smallest C with sum_j max(|w_j|-C,0) <= eps*sum|w|.
class ExcessUnweighter : public SortedWeightUnweighter, RegistrableUnweighter<ExcessUnweighter> {
  public:
    ExcessUnweighter(const YAML::Node &);
    static std::unique_ptr<Unweighter> Construct(const YAML::Node &);
    static std::string Name() { return "Excess"; }

  protected:
    double ComputeCap(const std::vector<double> &sorted_weights, double total) const override;
    std::string RuleTag() const override { return Name(); }
};

// Cap = smallest C with sum_{|w_j|>C} |w_j| <= eps*sum|w|.
class TailFractionUnweighter : public SortedWeightUnweighter,
                               RegistrableUnweighter<TailFractionUnweighter> {
  public:
    TailFractionUnweighter(const YAML::Node &);
    static std::unique_ptr<Unweighter> Construct(const YAML::Node &);
    static std::string Name() { return "TailFraction"; }

  protected:
    double ComputeCap(const std::vector<double> &sorted_weights, double total) const override;
    std::string RuleTag() const override { return Name(); }
};

} // namespace achilles

#endif
