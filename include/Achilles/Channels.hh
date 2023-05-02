#include "Achilles/Integrand.hh"
#include "Achilles/NuclearModel.hh"

// TODO: Turn this into a factory to reduce the number of includes
#include "Achilles/BeamMapper.hh"
#include "Achilles/FinalStateMapper.hh"
#include "Achilles/HadronicMapper.hh"
#include "Achilles/PhaseSpaceBuilder.hh"
#include "Achilles/PhaseSpaceMapper.hh"
#include "Achilles/QuasielasticTestMapper.hh"

#ifdef ACHILLES_SHERPA_INTERFACE
#include "plugins/Sherpa/Channels1.hh"
#include "plugins/Sherpa/Channels3.hh"
#endif

namespace achilles {

Channel<FourVector> BuildChannelTest(const YAML::Node &node, std::shared_ptr<Beam> beam);

template <typename T>
Channel<FourVector> BuildChannel(NuclearModel *model, size_t nlep, size_t nhad,
                                 std::shared_ptr<Beam> beam, const std::vector<double> &masses) {
    Channel<FourVector> channel;
    channel.mapping = PSBuilder(nlep, nhad)
                          .Beam(beam, masses)
                          .Hadron(model->PhaseSpace(), masses)
                          .FinalState(T::Name(), masses)
                          .build();
    AdaptiveMap map(channel.mapping->NDims(), 2);
    channel.integrator = Vegas(map, VegasParams{});
    return channel;
}

#ifdef ACHILLES_SHERPA_INTERFACE
Channel<FourVector> BuildGenChannel(NuclearModel *model, size_t nlep, size_t nhad,
                                    std::shared_ptr<Beam> beam,
                                    std::unique_ptr<PHASIC::Channels> final_state,
                                    const std::vector<double> &masses);
#endif

} // namespace achilles
