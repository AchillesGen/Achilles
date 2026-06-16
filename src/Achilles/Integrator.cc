#include "Achilles/Integrator.hh"
#include "Achilles/FourVector.hh"
#include "Achilles/Integrand.hh"
#include "Achilles/MultiChannel.hh"
#include "Achilles/Settings.hh"
#include "Achilles/Units.hh"

#include <algorithm>
#include <stdexcept>

#ifdef ACHILLES_ENABLE_PYTHON
#include "Achilles/FlowIntegrator.hh"
#endif

using achilles::Integrator;
using achilles::IntegratorType;

void achilles::Integrator::SampleForMaxWeight(Integrand<FourVector> &integrand, size_t nsamples) {
    const size_t batch_size = BatchSize();
    for(size_t drawn = 0; drawn < nsamples; drawn += batch_size) {
        const size_t n = std::min(batch_size, nsamples - drawn);
        auto points = GeneratePoints(integrand, n);
        auto weights = GenerateWeights(integrand, points);
        for(size_t i = 0; i < points.size(); ++i)
            if(weights[i] != 0) integrand(points[i], weights[i]);
    }
}

std::string_view achilles::ToString(IntegratorType type) {
    switch(type) {
    case IntegratorType::MultiChannel:
        return "MultiChannel";
    case IntegratorType::Flow:
        return "Flow";
    }
    throw std::runtime_error("Integrator: invalid IntegratorType value");
}

achilles::IntegratorType achilles::IntegratorTypeFromString(std::string_view name) {
    if(name == "MultiChannel") return IntegratorType::MultiChannel;
    if(name == "Flow") return IntegratorType::Flow;
    throw std::runtime_error("Integrator: unknown Options/Integrator/Type '" + std::string(name) +
                             "' (expected 'MultiChannel' or 'Flow').");
}

std::unique_ptr<Integrator> achilles::MakeIntegrator(const Settings &config, size_t ndims,
                                                     size_t nchannels) {
    UnitsEnum display_unit = UnitsEnum::nb;
    if(config.Exists("Main/DisplayUnit"))
        display_unit = ToUnitEnum(config.GetAs<std::string>("Main/DisplayUnit"));

    IntegratorType type = IntegratorType::MultiChannel;
    if(config.Exists("Options/Integrator/Type"))
        type = config.GetAs<IntegratorType>("Options/Integrator/Type");

    std::unique_ptr<Integrator> integrator;
    switch(type) {
    case IntegratorType::MultiChannel: {
        MultiChannelParams params;
        if(config.Exists("Options/Integrator/MultiChannel/Parameters"))
            params = config.GetAs<MultiChannelParams>("Options/Integrator/MultiChannel/Parameters");

        auto multichannel = std::make_unique<MultiChannel>(ndims, nchannels, params, display_unit);
        if(config.Exists("Options/Integrator/MultiChannel/Accuracy"))
            multichannel->SetRtol(config.GetAs<double>("Options/Integrator/MultiChannel/Accuracy"));
        integrator = std::move(multichannel);
        break;
    }
    case IntegratorType::Flow: {
#ifdef ACHILLES_ENABLE_PYTHON
        integrator = MakeFlowIntegrator(config, ndims, nchannels, display_unit);
        break;
#else
        throw std::runtime_error("Integrator: the Flow backend requires building Achilles with "
                                 "ACHILLES_ENABLE_PYTHON=ON.");
#endif
    }
    }

    if(!integrator) throw std::runtime_error("Integrator: unhandled integrator type.");

    // Common options shared by every backend.
    if(config.Exists("Options/Integrator/BatchSize"))
        integrator->SetBatchSize(config.GetAs<size_t>("Options/Integrator/BatchSize"));

    return integrator;
}
