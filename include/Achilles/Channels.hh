#include "Achilles/Integrand.hh"
#include "Achilles/NuclearModel.hh"
#include "Achilles/ProcessInfo.hh"

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
Channel<FourVector> BuildChannel(NuclearModel *model, const ProcessInfo &info,
                                 std::shared_ptr<Beam> beam, PID nuc_id,
                                 std::optional<double> gauge_boson_mass = std::nullopt) {
    Channel<FourVector> channel;
    channel.mapping = PSBuilder(info)
                          .Beam(beam)
                          .Hadron(model->PhaseSpace(nuc_id))
                          .FinalState(T::Name(), gauge_boson_mass)
                          .build();
    AdaptiveMap map(channel.mapping->NDims(), 100);
    channel.integrator = Vegas(map, VegasParams{});
    return channel;
}

#ifdef ACHILLES_SHERPA_INTERFACE
Channel<FourVector> BuildGenChannel(NuclearModel *model, const ProcessInfo &info,
                                    std::shared_ptr<Beam> beam,
                                    std::unique_ptr<PHASIC::Channels> final_state, PID nuc_id);
#endif

} // namespace achilles
