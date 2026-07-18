#include "Achilles/PyFlowSampler.hh"
#include "Achilles/FlowIntegrator.hh"
#include "Achilles/Settings.hh"

#include <mutex>
#include <set>

#include <pybind11/embed.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

namespace py = pybind11;

namespace {

// Ensure a Python interpreter exists process-wide.
void EnsureInterpreter() {
    static std::once_flag flag;
    static std::unique_ptr<py::scoped_interpreter> guard;
    std::call_once(flag, []() {
        if(!Py_IsInitialized()) guard = std::make_unique<py::scoped_interpreter>();
    });
}

// Convert a YAML scalar/sequence/map into the corresponding Python object so arbitrary
// run-card Flow parameters can be forwarded as keyword arguments. Scalars are probed as
// int, then float, then bool, then string.
py::object YamlToPy(const YAML::Node &node) {
    switch(node.Type()) {
    case YAML::NodeType::Scalar: {
        try {
            return py::int_(node.as<long long>());
        } catch(const YAML::Exception &) {}
        try {
            return py::float_(node.as<double>());
        } catch(const YAML::Exception &) {}
        try {
            return py::bool_(node.as<bool>());
        } catch(const YAML::Exception &) {}
        return py::str(node.as<std::string>());
    }
    case YAML::NodeType::Sequence: {
        py::list list;
        for(const auto &elem : node) list.append(YamlToPy(elem));
        return list;
    }
    case YAML::NodeType::Map: {
        py::dict dict;
        for(const auto &kv : node) dict[py::str(kv.first.as<std::string>())] = YamlToPy(kv.second);
        return dict;
    }
    default:
        return py::none();
    }
}

// Flatten an (n x ndims) vector-of-vectors into a 2D numpy array.
py::array_t<double> ToArray2D(const std::vector<std::vector<double>> &rows, size_t ndims) {
    const auto n = static_cast<py::ssize_t>(rows.size());
    py::array_t<double> arr({n, static_cast<py::ssize_t>(ndims)});
    auto view = arr.mutable_unchecked<2>();
    for(py::ssize_t i = 0; i < n; ++i)
        for(size_t j = 0; j < ndims; ++j)
            view(i, static_cast<py::ssize_t>(j)) = rows[static_cast<size_t>(i)][j];
    return arr;
}

py::array_t<long long> ToArray1D(const std::vector<size_t> &v) {
    py::array_t<long long> arr(static_cast<py::ssize_t>(v.size()));
    auto view = arr.mutable_unchecked<1>();
    for(size_t i = 0; i < v.size(); ++i)
        view(static_cast<py::ssize_t>(i)) = static_cast<long long>(v[i]);
    return arr;
}

py::array_t<double> ToArray1D(const std::vector<double> &v) {
    py::array_t<double> arr(static_cast<py::ssize_t>(v.size()));
    auto view = arr.mutable_unchecked<1>();
    for(size_t i = 0; i < v.size(); ++i) view(static_cast<py::ssize_t>(i)) = v[i];
    return arr;
}

std::vector<double>
ToVector(const py::array_t<double, py::array::c_style | py::array::forcecast> &arr) {
    auto view = arr.unchecked<1>();
    std::vector<double> out(static_cast<size_t>(view.shape(0)));
    for(py::ssize_t i = 0; i < view.shape(0); ++i) out[static_cast<size_t>(i)] = view(i);
    return out;
}

} // namespace

using achilles::PyFlowSampler;

PyFlowSampler::PyFlowSampler(const std::string &file, const std::string &class_name, size_t ndims,
                             size_t nchannels, py::dict params)
    : m_ndims{ndims} {
    EnsureInterpreter();
    py::gil_scoped_acquire gil;

    // load_flow loads the class (from the user file or the built-in module), instantiates
    // it, and validates that it implements the required protocol.
    py::object file_arg = py::none();
    if(!file.empty()) file_arg = py::cast(file);
    py::object loader = py::module_::import("Achilles.flow").attr("load_flow");
    m_flow = loader(file_arg, class_name, ndims, nchannels, **params);

    if(py::hasattr(m_flow, "requires_training"))
        m_requires_training = m_flow.attr("requires_training").cast<bool>();
}

achilles::FlowSampleBatch PyFlowSampler::Sample(size_t n) {
    py::gil_scoped_acquire gil;
    py::tuple result = m_flow.attr("sample")(n).cast<py::tuple>();

    auto channels =
        result[0].cast<py::array_t<long long, py::array::c_style | py::array::forcecast>>();
    auto rans = result[1].cast<py::array_t<double, py::array::c_style | py::array::forcecast>>();
    auto density = result[2].cast<py::array_t<double, py::array::c_style | py::array::forcecast>>();

    auto ch_view = channels.unchecked<1>();
    auto rn_view = rans.unchecked<2>();

    FlowSampleBatch batch;
    batch.channels.resize(n);
    batch.rans.assign(n, std::vector<double>(m_ndims));
    for(size_t i = 0; i < n; ++i) {
        batch.channels[i] = static_cast<size_t>(ch_view(static_cast<py::ssize_t>(i)));
        for(size_t j = 0; j < m_ndims; ++j)
            batch.rans[i][j] = rn_view(static_cast<py::ssize_t>(i), static_cast<py::ssize_t>(j));
    }
    batch.density = ToVector(density);
    return batch;
}

std::vector<double> PyFlowSampler::Density(const std::vector<size_t> &channels,
                                           const std::vector<std::vector<double>> &rans) {
    py::gil_scoped_acquire gil;
    py::object result = m_flow.attr("density")(ToArray1D(channels), ToArray2D(rans, m_ndims));
    return ToVector(result.cast<py::array_t<double, py::array::c_style | py::array::forcecast>>());
}

double PyFlowSampler::Train(const std::vector<size_t> &channels,
                            const std::vector<std::vector<double>> &rans,
                            const std::vector<double> &values) {
    py::gil_scoped_acquire gil;
    py::object loss =
        m_flow.attr("train")(ToArray1D(channels), ToArray2D(rans, m_ndims), ToArray1D(values));
    return loss.cast<double>();
}

void PyFlowSampler::Plot() const {
    py::gil_scoped_acquire gil;
    // Optional: only backends that record a training history expose this.
    if(py::hasattr(m_flow, "plot_loss")) m_flow.attr("plot_loss")();
}

void PyFlowSampler::Save(const std::string &path) const {
    py::gil_scoped_acquire gil;
    m_flow.attr("save")(path);
}

bool PyFlowSampler::Load(const std::string &path) {
    py::gil_scoped_acquire gil;
    try {
        return m_flow.attr("load")(path).cast<bool>();
    } catch(const py::error_already_set &) { return false; }
}

std::unique_ptr<achilles::Integrator> achilles::MakeFlowIntegrator(const Settings &config,
                                                                   size_t ndims, size_t nchannels,
                                                                   UnitsEnum unit) {
    EnsureInterpreter();
    py::gil_scoped_acquire gil;

    std::string file;
    if(config.Exists("Options/Integrator/Flow/File"))
        file = config.GetAs<std::string>("Options/Integrator/Flow/File");

    std::string class_name = "NormFlow";
    if(config.Exists("Options/Integrator/Flow/Class"))
        class_name = config.GetAs<std::string>("Options/Integrator/Flow/Class");

    // Every Flow key the integrator does not consume itself is forwarded to the flow
    // constructor as a keyword argument (architecture knobs, model/weights paths, ...).
    static const std::set<std::string> reserved = {
        "File", "Class", "Epochs", "NIterations", "NCalls", "Accuracy", "WeightsPath"};
    py::dict params;
    if(config.Exists("Options/Integrator/Flow")) {
        const YAML::Node flow_node = config["Options/Integrator/Flow"];
        for(const auto &kv : flow_node) {
            const auto key = kv.first.as<std::string>();
            if(reserved.count(key)) continue;
            params[py::str(key)] = YamlToPy(kv.second);
        }
    }

    auto sampler = std::make_unique<PyFlowSampler>(file, class_name, ndims, nchannels, params);
    // Tag the cache by the flow class so different flow backends for the same physics do
    // not share (and clobber) a cache entry.
    auto integrator = std::make_unique<FlowIntegrator>(ndims, nchannels, std::move(sampler), unit,
                                                       "Flow-" + class_name);

    if(config.Exists("Options/Integrator/Flow/NCalls"))
        integrator->SetNCalls(config.GetAs<size_t>("Options/Integrator/Flow/NCalls"));
    if(config.Exists("Options/Integrator/Flow/NIterations"))
        integrator->SetNIterations(config.GetAs<size_t>("Options/Integrator/Flow/NIterations"));
    if(config.Exists("Options/Integrator/Flow/Epochs"))
        integrator->SetMaxIterations(config.GetAs<size_t>("Options/Integrator/Flow/Epochs"));
    if(config.Exists("Options/Integrator/Flow/Accuracy"))
        integrator->SetRtol(config.GetAs<double>("Options/Integrator/Flow/Accuracy"));
    if(config.Exists("Options/Integrator/Flow/WeightsPath"))
        integrator->SetWeightsPath(
            config.GetAs<std::string>("Options/Integrator/Flow/WeightsPath"));

    return integrator;
}
