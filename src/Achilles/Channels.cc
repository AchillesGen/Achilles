#include "Achilles/Channels.hh"

using namespace achilles;

Channel<std::vector<FourVector>> achilles::BuildChannelTest(const YAML::Node &node, std::shared_ptr<Beam> beam) {
    Channel<std::vector<FourVector>> channel;
    channel.mapping = std::make_unique<QuasielasticTestMapper>(node, beam);
    AdaptiveMap map(channel.mapping->NDims(), 2);
    channel.integrator = Vegas(map, {});
    return channel;
}

#ifdef ACHILLES_SHERPA_INTERFACE
Channel<Event> achilles::BuildGenChannel(NuclearModel *model, const ProcessInfo &info,
                                              std::shared_ptr<Beam> beam,
                                              std::unique_ptr<PHASIC::Channels> final_state,
                                              PID nuc_id) {
    Channel<Event> channel;
    channel.mapping = PSBuilder(info)
                          .Beam(beam)
                          .Hadron(model->PhaseSpace(nuc_id))
                          .GenFinalState(std::move(final_state))
                          .build();
    AdaptiveMap map(channel.mapping->NDims(), 2);
    channel.integrator = Vegas(map, VegasParams{});
    return channel;
}
#endif
