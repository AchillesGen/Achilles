#ifndef ACHILLES_INTEGRATOR_HH
#define ACHILLES_INTEGRATOR_HH

#include <iosfwd>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "Achilles/Statistics.hh"

namespace achilles {

class FourVector;
class Settings;
template <typename T> class Integrand;

// Tag written as the first token of a serialized integrator stream so the matching
// backend can verify it is reconstructing the correct cache.
enum class IntegratorType { MultiChannel, Flow };

std::string_view ToString(IntegratorType);
IntegratorType IntegratorTypeFromString(std::string_view);

// Abstract integrator interface used by event generation.
class Integrator {
  public:
    Integrator() = default;
    Integrator(const Integrator &) = default;
    Integrator(Integrator &&) = default;
    Integrator &operator=(const Integrator &) = default;
    Integrator &operator=(Integrator &&) = default;
    virtual ~Integrator() = default;

    // One optimization/training iteration and the full optimization loop.
    virtual void operator()(Integrand<FourVector> &) = 0;
    virtual void Optimize(Integrand<FourVector> &) = 0;

    // Batched sampling. GeneratePoints returns n phase-space points; GenerateWeights
    // returns the phase-space weight for each.
    virtual std::vector<std::vector<FourVector>> GeneratePoints(Integrand<FourVector> &,
                                                                size_t n) = 0;
    virtual std::vector<double>
    GenerateWeights(Integrand<FourVector> &,
                    const std::vector<std::vector<FourVector>> &points) = 0;

    // Draws nsamples points from the current (fixed) sampling distribution, in batches of
    // BatchSize(), and evaluates the integrand on each nonzero-weight point.
    void SampleForMaxWeight(Integrand<FourVector> &integrand, size_t nsamples);

    // Number of points drawn per batched sampling call.
    size_t BatchSize() const { return m_batch_size; }
    void SetBatchSize(size_t n) { m_batch_size = std::max<size_t>(1, n); }

    // Status / control used by ProcessGroup.
    virtual bool HasResults() const = 0;
    virtual StatsData LastResult() const = 0;
    virtual bool NeedsOptimization(double rel_err) const = 0;
    virtual bool RequiresTraining() const { return true; }

    // Identifier for caching
    virtual std::string CacheTag() const = 0;
    virtual void UpdateSummary(bool) = 0;
    virtual void SetRtol(double) = 0;
    virtual void SetNCalls(size_t) = 0;

    virtual void SaveState(std::ostream &) const = 0;
    virtual bool LoadState(std::istream &) = 0;

  protected:
    static constexpr size_t batch_size_default = 1024;
    size_t m_batch_size{batch_size_default};
};

std::unique_ptr<Integrator> MakeIntegrator(const Settings &config, size_t ndims, size_t nchannels);

} // namespace achilles

namespace YAML {

template <> struct convert<achilles::IntegratorType> {
    static Node encode(const achilles::IntegratorType &rhs) {
        return Node(std::string(achilles::ToString(rhs)));
    }

    static bool decode(const Node &node, achilles::IntegratorType &rhs) {
        if(!node.IsScalar()) return false;
        rhs = achilles::IntegratorTypeFromString(node.as<std::string>());
        return true;
    }
};

} // namespace YAML

#endif
