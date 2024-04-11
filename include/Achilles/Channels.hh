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
Channel<FourVector> BuildChannel(NuclearModel *model, size_t nlep, size_t nhad, size_t nspec,
                                 std::shared_ptr<Beam> beam, const std::vector<double> &masses,
                                 PID nuc_id,
                                 std::optional<double> gauge_boson_mass = std::nullopt) {
    Channel<FourVector> channel;
    channel.mapping = PSBuilder(nlep, nhad, nspec)
                          .Beam(beam, masses)
                          .Hadron(model->PhaseSpace(nuc_id), masses)
                          .FinalState(T::Name(), masses, gauge_boson_mass)
                          .build();
    AdaptiveMap map(channel.mapping->NDims(), 2);
    channel.integrator = Vegas(map, VegasParams{});
    return channel;
}

#ifdef ACHILLES_SHERPA_INTERFACE
Channel<FourVector> BuildGenChannel(NuclearModel *model, size_t nlep, size_t nhad, size_t nspec,
                                    std::shared_ptr<Beam> beam,
                                    std::unique_ptr<PHASIC::Channels> final_state,
                                    const std::vector<double> &masses, PID nuc_id);
#endif

} // namespace achilles
