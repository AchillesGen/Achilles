#ifndef HARD_SCATTERING_ENUM_HH
#define HARD_SCATTERING_ENUM_HH

namespace achilles {

enum class HardScatteringType {
    None = -1,
    Coherent,
    Quasielastic,
    MesonExchangeCurrent,
    Interference_QE_MEC,
    Resonance,
    ShallowInelastic,
    DeepInelastic
};

}

#endif
