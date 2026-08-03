#include "Achilles/HadronicMapper.hh"
#include "Achilles/Constants.hh"
#include "Achilles/FourVector.hh"
#include "Achilles/ThreeVector.hh"
#include "Achilles/Event.hh"
#include "spdlog/spdlog.h"

using achilles::CoherentMapper;
using achilles::IntfSpectralMapper;
using achilles::QESpectralMapper;
using achilles::DISSingleNucleonMapper;
using achilles::DISNucleusMapper;

CoherentMapper::CoherentMapper(const ProcessInfo &info)
    : HadronicBeamMapper(info) {
    m_mass = ParticleInfo(info.m_hadronic.first[0]).Mass();
}

void CoherentMapper::GeneratePoint(Event& event, const std::vector<double> &) {
    event.addHadronIn({m_mass,0,0,0},ParticleStatus::target);
    Mapper<Event>::Print(__PRETTY_FUNCTION__, event, {});
}

double CoherentMapper::GenerateWeight(const Event&, std::vector<double> &) {
    return 1;
}


QESpectralMapper::QESpectralMapper(const ProcessInfo &info)
    : HadronicBeamMapper(info) {
    SetMasses(info.Masses());
}

void QESpectralMapper::GeneratePoint(Event& event,
                                     const std::vector<double> &rans) {
	FourVector& lepIn=event.LeptonsIn()[0].Momentum();
	// Generate inital nucleon state
	double radical =
	    pow(lepIn.E(), 2) + 2 * lepIn.E() * Constant::mN + Constant::mN2 - Smin();
	if(radical < 0) radical = 0;
	double pmin = lepIn.E() - sqrt(radical);
	double pmax = lepIn.E() + sqrt(radical);
	pmin = pmin < 0 ? 0 : pmin;
	pmax = pmax > 800 ? 800 : pmax;
	double dp = pmax - pmin;
	const double mom = dp * rans[0] + pmin;
	double cosT_max = (2 * lepIn.E() * Constant::mN + Constant::mN2 - mom * mom - Smin()) /
	                  (2 * lepIn.E() * mom);
	cosT_max = cosT_max > 1 ? 1 : cosT_max < -1 ? -1 : cosT_max;
	const double cosT = (cosT_max + 1) * rans[1] - 1;
	const double sinT = sqrt(1 - cosT * cosT);
	const double phi = dPhi * rans[2];
	ThreeVector pmom = {mom * sinT * cos(phi), mom * sinT * sin(phi), mom * cosT};

	const double det = pow(lepIn.E(), 2) + mom * mom + 2 * pmom * lepIn.Vec3() + Smin();
	double emax = Constant::mN + lepIn.E() - sqrt(det);
	emax = std::min(emax, Constant::mN - mom);
	emax = emax > 400 ? 400 : emax;
	const double energy = emax * rans[3] - 1e-8;
	// if(emax < 0) energy = emax - 1;
	// const double energy = dE*rans[3];

	// double cosT_max = (Constant::mN2 + energy*energy - 2*Constant::mN*energy
	// - mom*mom + 2*point[1].E()*Constant::mN - 2*point[1].E()*energy -
	// Smin())/(2*mom*point[1].P()); cosT_max = cosT_max > 1 ? 1 : cosT_max;
	// const double cosT = (cosT_max + 1) * rans[1] - 1;
	// const double sinT = sqrt(1 - cosT*cosT);

	event.addHadronIn({Constant::mN-energy, mom*sinT*cos(phi), mom*sinT*sin(phi), mom*cosT});
#ifdef ACHILLES_EVENT_DETAILS
    Mapper<Event>::Print(__PRETTY_FUNCTION__, event, rans);
    spdlog::trace("  lepIn = {}", lepIn);
    spdlog::trace("  dp = {}", dp);
    spdlog::trace("  cosT_max = {}", cosT_max);
    spdlog::trace("  cosT = {}", cosT);
    spdlog::trace("  mom = {}", mom);
    spdlog::trace("  energy = {}", energy);
    spdlog::trace("  emax = {}", emax);
    spdlog::trace("  s = {}", (lepIn + hadIn).M2());
    spdlog::trace("  s_min = {}", Smin());
#endif
}

double QESpectralMapper::GenerateWeight(const Event& event,
                                        std::vector<double> &rans) {
	const FourVector& lepIn=event.LeptonsIn()[0].Momentum(),hadIn=event.HadronsIn()[0].Momentum();
    double radical = pow(lepIn.E(), 2) + 2 * lepIn.E() * Constant::mN + Constant::mN2 - Smin();
    if(radical < 0) radical = 0;
    double pmin = lepIn.E() - sqrt(radical);
    double pmax = lepIn.E() + sqrt(radical);
    pmin = pmin < 0 ? 0 : pmin;
    pmax = pmax > 800 ? 800 : pmax;
    double dp = pmax - pmin;
    rans[0] = (hadIn.P() - pmin) / dp;
    double cosT_max = (2 * lepIn.E() * Constant::mN + Constant::mN2 - hadIn.P2() - Smin()) /
                      (2 * lepIn.E() * hadIn.P());
    cosT_max = cosT_max > 1 ? 1 : cosT_max < -1 ? -1 : cosT_max;
    double dCos = (cosT_max + 1);
    rans[1] = (hadIn.CosTheta() + 1) / dCos;
    rans[2] = hadIn.Phi() / dPhi;

    const double det =
        pow(lepIn.E(), 2) + hadIn.P2() + 2 * hadIn.Vec3() * lepIn.Vec3() + Smin();
    double emax = Constant::mN + lepIn.E() - sqrt(det);
    emax = std::min(emax, Constant::mN - hadIn.P());
    emax = emax > 400 ? 400 : emax;
    const double energy = Constant::mN - hadIn.E();
    // if(energy < 0) return std::numeric_limits<double>::infinity();
    const double dE = emax;
    rans[3] = (energy + 1e-8) / emax;
    // rans[3] = (Constant::mN - hadIn.E())/dE;

    // double cosT_max = (hadIn.M2() +
    // 2*point[1].E()*hadIn.E() -
    // Smin())/(2*hadIn.P()*point[1].P()); cosT_max = cosT_max > 1
    // ? 1 : cosT_max; const double dCos = (cosT_max + 1);
    // rans[1] = (hadIn.CosTheta() + 1) / dCos;

    double wgt = 1.0 / hadIn.P2() / dp / dCos / dPhi / dE;
#ifdef ACHILLES_EVENT_DETAILS
    Mapper<Event>::Print(__PRETTY_FUNCTION__, event, rans);
    spdlog::trace("  Weight: {}", wgt);
    spdlog::trace("  dp: {}", dp);
    spdlog::trace("  dCos: {}", dCos);
    spdlog::trace("  dPhi: {}", dPhi);
    spdlog::trace("  dE: {}", dE);
#endif

    return wgt;
}


IntfSpectralMapper::IntfSpectralMapper(const ProcessInfo &info)
    : HadronicBeamMapper(info) {
    SetMasses(info.Masses());
}

void IntfSpectralMapper::GeneratePoint(Event& event,
                                       const std::vector<double> &rans) {
	FourVector& lepIn=event.LeptonsIn()[0].Momentum();
    // Generate spectator momentum from a flat distribution
    double p2 = dp2 * rans[4];
    double cosT2 = dCos2 * rans[5] - 1.;
    double sinT2 = sqrt(1 - cosT2 * cosT2);
    double phi2 = dPhi * rans[6];

    ThreeVector pmom2 = {p2 * sinT2 * cos(phi2), p2 * sinT2 * sin(phi2), p2 * cosT2};
    double E2 = sqrt(p2 * p2 + Constant::mN2);

    event.addSpectator({E2, pmom2[0], pmom2[1], pmom2[2]});

    // Generate inital nucleon state
    double pmin = lepIn.E() - sqrt(pow(lepIn.E(), 2) + 2 * lepIn.E() * Constant::mN +
                                      Constant::mN2 - Smin());
    double pmax = lepIn.E() + sqrt(pow(lepIn.E(), 2) + 2 * lepIn.E() * Constant::mN +
                                      Constant::mN2 - Smin());
    pmin = pmin < 0 ? 0 : pmin;
    pmax = pmax > 800 ? 800 : pmax;
    double dp = pmax - pmin;
    const double mom = dp * rans[0] + pmin;
    double cosT_max = (2 * lepIn.E() * Constant::mN + Constant::mN2 - mom * mom - Smin()) /
                      (2 * lepIn.E() * mom);
    cosT_max = cosT_max > 1 ? 1 : cosT_max;
    const double cosT = (cosT_max + 1) * rans[1] - 1;
    const double sinT = sqrt(1 - cosT * cosT);
    const double phi = dPhi * rans[2];
    ThreeVector pmom = {mom * sinT * cos(phi), mom * sinT * sin(phi), mom * cosT};

    const double det = pow(lepIn.E(), 2) + mom * mom + 2 * pmom * lepIn.Vec3() + Smin();
    double emax = Constant::mN + lepIn.E() - sqrt(det);
    emax = std::min(emax, Constant::mN - mom);
    emax = emax > 400 ? 400 : emax;
    const double energy = emax * rans[3] - 1e-8;
    // if(emax < 0) energy = emax - 1;
    // const double energy = dE*rans[3];

    // double cosT_max = (Constant::mN2 + energy*energy - 2*Constant::mN*energy
    // - mom*mom + 2*point[1].E()*Constant::mN - 2*point[1].E()*energy -
    // Smin())/(2*mom*point[1].P()); cosT_max = cosT_max > 1 ? 1 : cosT_max;
    // const double cosT = (cosT_max + 1) * rans[1] - 1;
    // const double sinT = sqrt(1 - cosT*cosT);

    event.addHadronIn({Constant::mN-energy, mom*sinT*cos(phi), mom*sinT*sin(phi), mom*cosT});
#ifdef ACHILLES_EVENT_DETAILS
    Mapper<Event>::Print(__PRETTY_FUNCTION__, event, rans);
    spdlog::trace("  lepIn = {}", lepIn);
    spdlog::trace("  dp = {}", dp);
    spdlog::trace("  cosT_min = {}", cosT_max);
    spdlog::trace("  cosT_max = {}", cosT_max);
    spdlog::trace("  cosT = {}", cosT);
    spdlog::trace("  mom = {}", mom);
    spdlog::trace("  energy = {}", energy);
    spdlog::trace("  emax = {}", emax);
    spdlog::trace("  s = {}", (lepIn + point[1]).M2());
    spdlog::trace("  s_min = {}", Smin());
#endif
}

double IntfSpectralMapper::GenerateWeight(const Event& event,
                                          std::vector<double> &rans) {
	const FourVector& lepIn=event.LeptonsIn()[0].Momentum(),hadIn=event.HadronsIn()[0].Momentum(),
					  spectator=event.Spectators()[0].Momentum();

    double pmin = lepIn.E() - sqrt(pow(lepIn.E(), 2) + 2 * lepIn.E() * Constant::mN +
                                      Constant::mN2 - Smin());
    double pmax = lepIn.E() + sqrt(pow(lepIn.E(), 2) + 2 * lepIn.E() * Constant::mN +
                                      Constant::mN2 - Smin());
    pmin = pmin < 0 ? 0 : pmin;
    pmax = pmax > 800 ? 800 : pmax;
    double dp = pmax - pmin;
    rans[0] = (hadIn.P() - pmin) / dp;
    double cosT_max = (2 * lepIn.E() * Constant::mN + Constant::mN2 - hadIn.P2() - Smin()) /
                      (2 * lepIn.E() * hadIn.P());
    cosT_max = cosT_max > 1 ? 1 : cosT_max;
    double dCos = (cosT_max + 1);
    rans[1] = (hadIn.CosTheta() + 1) / dCos;
    rans[2] = hadIn.Phi() / dPhi;

    const double det =
        pow(lepIn.E(), 2) + hadIn.P2() + 2 * hadIn.Vec3() * lepIn.Vec3() + Smin();
    double emax = Constant::mN + lepIn.E() - sqrt(det);
    emax = std::min(emax, Constant::mN - hadIn.P());
    emax = emax > 400 ? 400 : emax;
    const double energy = Constant::mN - hadIn.E();
    // if(energy < 0) return std::numeric_limits<double>::infinity();
    const double dE = emax;
    rans[3] = (energy + 1e-8) / emax;

    rans[4] = spectator.P() / dp2;
    rans[5] = (spectator.CosTheta() + 1) / dCos2;
    rans[6] = spectator.Phi() / dPhi;

    // rans[3] = (Constant::mN - hadIn.E())/dE;

    // double cosT_max = (hadIn.M2() +
    // 2*point[1].E()*hadIn.E() -
    // Smin())/(2*hadIn.P()*point[1].P()); cosT_max = cosT_max > 1
    // ? 1 : cosT_max; const double dCos = (cosT_max + 1);
    // rans[1] = (hadIn.CosTheta() + 1) / dCos;

    double wgt =
        1.0 / hadIn.P2() / dp / dCos / dPhi / dE / spectator.P2() / dp2 / dCos2 / dPhi;
#ifdef ACHILLES_EVENT_DETAILS
    Mapper<Event>::Print(__PRETTY_FUNCTION__, event, rans);
    spdlog::trace("  Weight: {}", wgt);
    spdlog::trace("  dp: {}", dp);
    spdlog::trace("  dCos: {}", dCos);
    spdlog::trace("  dPhi: {}", dPhi);
    spdlog::trace("  dE: {}", dE);
#endif

    return wgt;
}


static achilles::FourVector getQuarkMomentum(achilles::FourVector& pLepton,achilles::FourVector& pHadron,double x) {
	// ~~ Boost to COM Frame ~~
	achilles::ThreeVector com=(pLepton+pHadron).BoostVector();
	achilles::FourVector pHadBoosted=pHadron.Boost(com);

	// Rotated hadron momentum is always [E,0,0,pZ=|p|]
	// Instead of rotating I'll just call P() to act as pZ

	// ~~ Lightcone Coords ~~
	// p+_hadron = (Ehad + pZhad)/sqrt(2)
	// p+ = x * p+_hadron
	// pQuark=[(p+ + p-)/sqrt(2),0,0,(p+ - p-)/sqrt(2)] (regular coords)
	// p- = mQuark^2/2p+ = 0 for massless quark
	// => Equark = pZquark (rotated) = p+/sqrt(2) = x*(Ehad+pZhad)/2
	double Eq=x*(pHadBoosted.E()+pHadBoosted.P())/2.0;
	// => pQuark (COM frame, rotated) = [Eq,0,0,Eq]

	// ~~ Rotate/Boost Back ~~
	// pQuark (COM frame, unrotated) = [Eq,Eq*(original unit vector)]
	// Want to return in Lab Frame (so boost by -com)
	return achilles::FourVector(pHadBoosted.Vec3().Unit()*Eq,Eq).Boost(-com);
}

static double getX(const achilles::FourVector& pLepton,const achilles::FourVector& pHadron,const achilles::FourVector& pQuark) {
	// Reversing the previous function
	achilles::ThreeVector com=(pLepton+pHadron).BoostVector();
	achilles::FourVector pHadBoosted=pHadron.Boost(com);
	return 2.0*pQuark.Boost(com).E()/(pHadBoosted.E()+pHadBoosted.P());
}

static double xMin=1e-10,xMax=1.0;

DISSingleNucleonMapper::DISSingleNucleonMapper(const ProcessInfo& info): HadronicBeamMapper(info) {
	SetMasses(info.Masses());
	// Hadron is always stationary in this model
	pHadron={ParticleInfo(info.m_hadronic.first[0]).Mass(),0,0,0};
}

void DISSingleNucleonMapper::GeneratePoint(Event& event,const std::vector<double>& rans) {
	// Generate initial quark state
	double x=xMin*pow(xMax/xMin,rans[0]);
	event.addHadronIn(getQuarkMomentum(event.LeptonsIn()[0].Momentum(),pHadron,x));
}

double DISSingleNucleonMapper::GenerateWeight(const Event& event,std::vector<double>& rans) {
	const FourVector& hadIn=event.HadronsIn()[0].Momentum();
	double x=getX(event.LeptonsIn()[0].Momentum(),pHadron,hadIn);
	rans[0]=log(x/xMin)/log(xMax/xMin);
	return 1.0/(hadIn.P2()*(x*log(xMax/xMin)));
}


DISNucleusMapper::DISNucleusMapper(const ProcessInfo& info): HadronicBeamMapper(info) {
	SetMasses(info.Masses());
}

void DISNucleusMapper::GeneratePoint(Event& event,const std::vector<double>& rans) {
	FourVector& lepIn=event.LeptonsIn()[0].Momentum();
	// Generate inital nucleon state
	double radical =
		pow(lepIn.E(), 2) + 2 * lepIn.E() * Constant::mN + Constant::mN2 - Smin();
	if(radical < 0) radical = 0;
	double pmin = lepIn.E() - sqrt(radical);
	double pmax = lepIn.E() + sqrt(radical);
	pmin = pmin < 0 ? 0 : pmin;
	pmax = pmax > 800 ? 800 : pmax;
	double dp = pmax - pmin;
	const double mom = dp * rans[0] + pmin;
	double cosT_max = (2 * lepIn.E() * Constant::mN + Constant::mN2 - mom * mom - Smin()) /
					  (2 * lepIn.E() * mom);
	cosT_max = cosT_max > 1 ? 1 : cosT_max < -1 ? -1 : cosT_max;
	const double cosT = (cosT_max + 1) * rans[1] - 1;
	const double sinT = sqrt(1 - cosT * cosT);
	const double phi = dPhi * rans[2];
	ThreeVector pmom = {mom * sinT * cos(phi), mom * sinT * sin(phi), mom * cosT};

	const double det = pow(lepIn.E(), 2) + mom * mom + 2 * pmom * lepIn.Vec3() + Smin();
	double emax = Constant::mN + lepIn.E() - sqrt(det);
	emax = std::min(emax, Constant::mN - mom);
	emax = emax > 400 ? 400 : emax;
	const double energy = emax * rans[3] - 1e-8;

	// pHadron
	FourVector hadIn={Constant::mN - energy, mom * sinT * cos(phi), mom * sinT * sin(phi), mom * cosT};
	event.addHadronIn(hadIn);

	double x=xMin*pow(xMax/xMin,rans[4]);

	// pQuark
	event.addHadronIn(getQuarkMomentum(lepIn,hadIn,x));
}

double DISNucleusMapper::GenerateWeight(const Event& event,std::vector<double>& rans) {
	const FourVector& lepIn=event.LeptonsIn()[0].Momentum(),
					  hadIn=event.HadronsIn()[0].Momentum(),
					  quarkIn=event.HadronsIn()[1].Momentum();
	double x=getX(lepIn,hadIn,quarkIn);

	double radical =
	    pow(lepIn.E(), 2) + 2 * lepIn.E() * Constant::mN + Constant::mN2 - Smin();
	if(radical < 0) radical = 0;
	double pmin = lepIn.E() - sqrt(radical);
	double pmax = lepIn.E() + sqrt(radical);
	pmin = pmin < 0 ? 0 : pmin;
	pmax = pmax > 800 ? 800 : pmax;
	double dp = pmax - pmin;
	rans[0] = (hadIn.P() - pmin) / dp;
	double cosT_max = (2 * lepIn.E() * Constant::mN + Constant::mN2 - hadIn.P2() - Smin()) /
	                  (2 * lepIn.E() * hadIn.P());
	cosT_max = cosT_max > 1 ? 1 : cosT_max < -1 ? -1 : cosT_max;
	double dCos = (cosT_max + 1);
	rans[1] = (hadIn.CosTheta() + 1) / dCos;
	rans[2] = hadIn.Phi() / dPhi;

	const double det =
	    pow(lepIn.E(), 2) + hadIn.P2() + 2 * hadIn.Vec3() * lepIn.Vec3() + Smin();
	double emax = Constant::mN + lepIn.E() - sqrt(det);
	emax = std::min(emax, Constant::mN - hadIn.P());
	emax = emax > 400 ? 400 : emax;
	const double energy = Constant::mN - hadIn.E();
	const double dE = emax;
	rans[3] = (energy + 1e-8) / emax;

	rans[4]=log(x/xMin)/log(xMax/xMin);

    return 1.0/(quarkIn.P2()*dp*dCos*dPhi*dE*(x*log(xMax/xMin)));
}