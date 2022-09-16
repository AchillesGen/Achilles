#include "Achilles/Integrand.hh"
#include "Achilles/NuclearModel.hh"

// TODO: Turn this into a factory to reduce the number of includes
#include "Achilles/PhaseSpaceBuilder.hh"
#include "Achilles/BeamMapper.hh"
#include "Achilles/HadronicMapper.hh"
#include "Achilles/FinalStateMapper.hh"
#include "Achilles/PhaseSpaceMapper.hh"
#include "Achilles/QuasielasticTestMapper.hh"

#ifdef ENABLE_BSM
#include "plugins/Sherpa/Channels1.hh"
#include "plugins/Sherpa/Channels3.hh"
#include "plugins/Sherpa/SherpaMEs.hh"
#endif

namespace achilles {

Channel<FourVector> BuildChannelTest(const YAML::Node &node, std::shared_ptr<Beam> beam);

template<typename T>
Channel<FourVector> BuildChannel(NuclearModel *model, size_t nlep, size_t nhad,
                                 std::shared_ptr<Beam> beam,
                                 const std::vector<double> &masses) {
    Channel<FourVector> channel;
    channel.mapping = PSBuilder(nlep, nhad).Beam(beam, 1)
                                           .Hadron(model -> PhaseSpace(), masses)
                                           .FinalState(T::Name(), masses).build();
    AdaptiveMap map(channel.mapping -> NDims(), 2);
    channel.integrator = Vegas(map, VegasParams{});
    return channel;
}

#ifdef ENABLE_BSM
Channel<FourVector> BuildGenChannel(NuclearModel *model, size_t nlep, size_t nhad,
                                    std::shared_ptr<Beam> beam,
                                    std::unique_ptr<PHASIC::Channels> final_state,
                                    const std::vector<double> &masses);
#endif

}
