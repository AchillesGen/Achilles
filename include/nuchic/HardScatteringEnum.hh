#ifndef HARD_SCATTERING_ENUM_HH
#define HARD_SCATTERING_ENUM_HH

namespace nuchic {

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
