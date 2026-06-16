#ifndef ACHILLES_FLOWSAMPLER_HH
#define ACHILLES_FLOWSAMPLER_HH

#include <cstddef>
#include <string>
#include <vector>

namespace achilles {

struct FlowSampleBatch {
    std::vector<size_t> channels;
    std::vector<std::vector<double>> rans;
    std::vector<double> density;
};

// Abstract sampling backend behind FlowIntegrator. Decouples the integrator's
// orchestration from any particular flow implementation (and from Python/torch), so the
// integrator can be unit-tested.
class FlowSampler {
  public:
    FlowSampler() = default;
    FlowSampler(const FlowSampler &) = default;
    FlowSampler(FlowSampler &&) = default;
    FlowSampler &operator=(const FlowSampler &) = default;
    FlowSampler &operator=(FlowSampler &&) = default;
    virtual ~FlowSampler() = default;

    // Draw n samples (channel, rans, density) from the current flow.
    virtual FlowSampleBatch Sample(size_t n) = 0;

    // Sampling density g(channel, rans) for externally supplied points. Used to build the
    // multichannel weight: the same point is queried under every channel.
    virtual std::vector<double> Density(const std::vector<size_t> &channels,
                                        const std::vector<std::vector<double>> &rans) = 0;

    // One optimization step on a batch of (channel, rans) labelled by integrand values.
    // Returns the loss for diagnostics.
    virtual double Train(const std::vector<size_t> &channels,
                         const std::vector<std::vector<double>> &rans,
                         const std::vector<double> &values) = 0;

    // Whether this flow still needs training.
    virtual bool RequiresTraining() const { return true; }

    // Persist / restore the learned parameters (e.g. torch weights) to a path.
    virtual void Save(const std::string &path) const = 0;
    virtual bool Load(const std::string &path) = 0;
};

} // namespace achilles

#endif
