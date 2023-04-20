#include "Achilles/Channels.hh"
#include "Achilles/Particle.hh"

using namespace achilles;

Channel<FourVector> achilles::BuildChannelTest(const YAML::Node &node, std::shared_ptr<Beam> beam) {
    Channel<FourVector> channel;
    channel.mapping = std::make_unique<QuasielasticTestMapper>(node, beam);
    AdaptiveMap map(channel.mapping -> NDims(), 2);
    channel.integrator = Vegas(map, {});
    return channel;
}

#ifdef ENABLE_BSM
Channel<FourVector> achilles::BuildGenChannel(NuclearModel *model, size_t nlep, size_t nhad,
                                              std::shared_ptr<Beam> beam,
                                              std::unique_ptr<PHASIC::Channels> final_state,
                                              const std::vector<double> &masses) {
    Channel<FourVector> channel;
    channel.mapping = PSBuilder(nlep, nhad).Beam(beam, 1, masses)
                                           .Hadron(model -> PhaseSpace(), masses)
                                           .GenFinalState(std::move(final_state)).build();
    AdaptiveMap map(channel.mapping -> NDims(), 2);
    channel.integrator = Vegas(map, VegasParams{});
    return channel;
}
#endif
