#ifndef ACHILLES_PYFLOWSAMPLER_HH
#define ACHILLES_PYFLOWSAMPLER_HH

#include <string>

#include <pybind11/pybind11.h>

#include "Achilles/FlowSampler.hh"

namespace achilles {

// FlowSampler backed by a Python flow object. The object may be one of the built-in
// classes in python/Achilles/flow, or any class the user supplies in their own .py file
// (selected by File + Class in the run card). Arrays are exchanged as numpy.
class PyFlowSampler : public FlowSampler {
  public:
    // `file` is an optional path to a user module (empty => use the built-in
    // Achilles.flow classes). `params` are forwarded as keyword arguments to the flow
    // class constructor.
    PyFlowSampler(const std::string &file, const std::string &class_name, size_t ndims,
                  size_t nchannels, pybind11::dict params);

    FlowSampleBatch Sample(size_t n) override;
    std::vector<double> Density(const std::vector<size_t> &channels,
                                const std::vector<std::vector<double>> &rans) override;
    double Train(const std::vector<size_t> &channels, const std::vector<std::vector<double>> &rans,
                 const std::vector<double> &values) override;
    void Plot() const override;
    bool RequiresTraining() const override { return m_requires_training; }
    void Save(const std::string &path) const override;
    bool Load(const std::string &path) override;

  private:
    size_t m_ndims;
    bool m_requires_training{true};
    pybind11::object m_flow;
};

} // namespace achilles

#endif
