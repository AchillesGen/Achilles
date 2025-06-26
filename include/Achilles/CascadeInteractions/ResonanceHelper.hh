#ifndef RESONANCEHELPER_HH
#define RESONANCEHELPER_HH

#include <cstddef>

namespace achilles {

class PID;
class Particle;
class Random;

namespace resonance {
// NN to NDelta and NDelta to NN
double Pcm2(double, double, double);
double DSigmaDM(bool iresonance, double sqrts, double mdelta, PID delta_id);
double MatNN2NDelta(double t, double u, double mdelta, PID delta_id);

// NDelta to NDelta
double DSigmaND2ND(double sqrts, double mn1, double mn2, double mu1, double mu2, double spectral);
double MatNDelta2NDelta(double t, double mu1, double mu2);

// Resonance mass and widths
double BreitWignerSpectral(PID id, double mu2);
double GetEffectiveWidth(PID id, double mass, double mass1, double mass2, size_t angular_mom);
double GenerateMass(const Particle &p1, const Particle &p2, PID res, PID other, Random &ran,
                    double smax);

} // namespace resonance

} // namespace achilles

#endif // RESONANCEHELPER_HH
