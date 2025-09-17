#ifndef ACHILLES_CLEBSCHGORDAN
#define ACHILLES_CLEBSCHGORDAN

namespace achilles {

struct SpinState {
    int total, zaxis;
    SpinState() = default;
    SpinState(double _total, double _zaxis)
        : total{static_cast<int>(2 * _total)}, zaxis{static_cast<int>(2 * _zaxis)} {}
};

double ClebschGordan(SpinState, SpinState, SpinState);

} // namespace achilles

#endif
