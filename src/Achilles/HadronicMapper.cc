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

CoherentMapper::CoherentMapper(const ProcessInfo &info, size_t idx)
    : HadronicBeamMapper(info, idx) {
    m_mass = ParticleInfo(info.m_hadronic.first[0]).Mass();
}

void CoherentMapper::GeneratePoint(Event& event, const std::vector<double> &) {
    event.Momentum()[HadronIdx()] = {m_mass, 0, 0, 0};
    Mapper<Event>::Print(__PRETTY_FUNCTION__, event, {});
}

double CoherentMapper::GenerateWeight(const Event&, std::vector<double> &) {
    return 1;
}


QESpectralMapper::QESpectralMapper(const ProcessInfo &info, size_t idx)
    : HadronicBeamMapper(info, idx) {
    SetMasses(info.Masses());
}

void QESpectralMapper::GeneratePoint(Event& event,
                                     const std::vector<double> &rans) {
	std::vector<FourVector>& point=event.Momentum();
    // Generate inital nucleon state
    double radical =
        pow(point[0].E(), 2) + 2 * point[0].E() * Constant::mN + Constant::mN2 - Smin();
    if(radical < 0) radical = 0;
    double pmin = point[0].E() - sqrt(radical);
    double pmax = point[0].E() + sqrt(radical);
    pmin = pmin < 0 ? 0 : pmin;
    pmax = pmax > 800 ? 800 : pmax;
    double dp = pmax - pmin;
    const double mom = dp * rans[0] + pmin;
    double cosT_max = (2 * point[0].E() * Constant::mN + Constant::mN2 - mom * mom - Smin()) /
                      (2 * point[0].E() * mom);
    cosT_max = cosT_max > 1 ? 1 : cosT_max < -1 ? -1 : cosT_max;
    const double cosT = (cosT_max + 1) * rans[1] - 1;
    const double sinT = sqrt(1 - cosT * cosT);
    const double phi = dPhi * rans[2];
    ThreeVector pmom = {mom * sinT * cos(phi), mom * sinT * sin(phi), mom * cosT};

    const double det = pow(point[0].E(), 2) + mom * mom + 2 * pmom * point[0].Vec3() + Smin();
    double emax = Constant::mN + point[0].E() - sqrt(det);
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

    point[HadronIdx()] = {Constant::mN - energy, mom * sinT * cos(phi), mom * sinT * sin(phi),
                          mom * cosT};
#ifdef ACHILLES_EVENT_DETAILS
    Mapper<Event>::Print(__PRETTY_FUNCTION__, event, rans);
    spdlog::trace("  point[0] = {}", point[0]);
    spdlog::trace("  dp = {}", dp);
    spdlog::trace("  cosT_max = {}", cosT_max);
    spdlog::trace("  cosT = {}", cosT);
    spdlog::trace("  mom = {}", mom);
    spdlog::trace("  energy = {}", energy);
    spdlog::trace("  emax = {}", emax);
    spdlog::trace("  s = {}", (point[0] + point[HadronIdx()]).M2());
    spdlog::trace("  s_min = {}", Smin());
#endif
}

double QESpectralMapper::GenerateWeight(const Event& event,
                                        std::vector<double> &rans) {
	const std::vector<FourVector>& point=event.Momentum();
    double radical =
        pow(point[0].E(), 2) + 2 * point[0].E() * Constant::mN + Constant::mN2 - Smin();
    if(radical < 0) radical = 0;
    double pmin = point[0].E() - sqrt(radical);
    double pmax = point[0].E() + sqrt(radical);
    pmin = pmin < 0 ? 0 : pmin;
    pmax = pmax > 800 ? 800 : pmax;
    double dp = pmax - pmin;
    rans[0] = (point[HadronIdx()].P() - pmin) / dp;
    double cosT_max = (2 * point[0].E() * Constant::mN + Constant::mN2 - point[HadronIdx()].P2() - Smin()) /
                      (2 * point[0].E() * point[HadronIdx()].P());
    cosT_max = cosT_max > 1 ? 1 : cosT_max < -1 ? -1 : cosT_max;
    double dCos = (cosT_max + 1);
    rans[1] = (point[HadronIdx()].CosTheta() + 1) / dCos;
    rans[2] = point[HadronIdx()].Phi() / dPhi;

    const double det =
        pow(point[0].E(), 2) + point[HadronIdx()].P2() + 2 * point[HadronIdx()].Vec3() * point[0].Vec3() + Smin();
    double emax = Constant::mN + point[0].E() - sqrt(det);
    emax = std::min(emax, Constant::mN - point[HadronIdx()].P());
    emax = emax > 400 ? 400 : emax;
    const double energy = Constant::mN - point[HadronIdx()].E();
    // if(energy < 0) return std::numeric_limits<double>::infinity();
    const double dE = emax;
    rans[3] = (energy + 1e-8) / emax;
    // rans[3] = (Constant::mN - point[HadronIdx()].E())/dE;

    // double cosT_max = (point[HadronIdx()].M2() +
    // 2*point[1].E()*point[HadronIdx()].E() -
    // Smin())/(2*point[HadronIdx()].P()*point[1].P()); cosT_max = cosT_max > 1
    // ? 1 : cosT_max; const double dCos = (cosT_max + 1);
    // rans[1] = (point[HadronIdx()].CosTheta() + 1) / dCos;

    double wgt = 1.0 / point[HadronIdx()].P2() / dp / dCos / dPhi / dE;
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


IntfSpectralMapper::IntfSpectralMapper(const ProcessInfo &info, size_t idx)
    : HadronicBeamMapper(info, idx) {
    SetMasses(info.Masses());
}

void IntfSpectralMapper::GeneratePoint(Event& event,
                                       const std::vector<double> &rans) {
	std::vector<FourVector>& point=event.Momentum();
    // Generate spectator momentum from a flat distribution
    double p2 = dp2 * rans[4];
    double cosT2 = dCos2 * rans[5] - 1.;
    double sinT2 = sqrt(1 - cosT2 * cosT2);
    double phi2 = dPhi * rans[6];

    ThreeVector pmom2 = {p2 * sinT2 * cos(phi2), p2 * sinT2 * sin(phi2), p2 * cosT2};
    double E2 = sqrt(p2 * p2 + Constant::mN2);

    point.back() = {E2, pmom2[0], pmom2[1], pmom2[2]};

    // Generate inital nucleon state
    double pmin = point[0].E() - sqrt(pow(point[0].E(), 2) + 2 * point[0].E() * Constant::mN +
                                      Constant::mN2 - Smin());
    double pmax = point[0].E() + sqrt(pow(point[0].E(), 2) + 2 * point[0].E() * Constant::mN +
                                      Constant::mN2 - Smin());
    pmin = pmin < 0 ? 0 : pmin;
    pmax = pmax > 800 ? 800 : pmax;
    double dp = pmax - pmin;
    const double mom = dp * rans[0] + pmin;
    double cosT_max = (2 * point[0].E() * Constant::mN + Constant::mN2 - mom * mom - Smin()) /
                      (2 * point[0].E() * mom);
    cosT_max = cosT_max > 1 ? 1 : cosT_max;
    const double cosT = (cosT_max + 1) * rans[1] - 1;
    const double sinT = sqrt(1 - cosT * cosT);
    const double phi = dPhi * rans[2];
    ThreeVector pmom = {mom * sinT * cos(phi), mom * sinT * sin(phi), mom * cosT};

    const double det = pow(point[0].E(), 2) + mom * mom + 2 * pmom * point[0].Vec3() + Smin();
    double emax = Constant::mN + point[0].E() - sqrt(det);
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

    point[HadronIdx()] = {Constant::mN - energy, mom * sinT * cos(phi), mom * sinT * sin(phi),
                          mom * cosT};
#ifdef ACHILLES_EVENT_DETAILS
    Mapper<Event>::Print(__PRETTY_FUNCTION__, event, rans);
    spdlog::trace("  point[0] = {}", point[0]);
    spdlog::trace("  dp = {}", dp);
    spdlog::trace("  cosT_min = {}", cosT_max);
    spdlog::trace("  cosT_max = {}", cosT_max);
    spdlog::trace("  cosT = {}", cosT);
    spdlog::trace("  mom = {}", mom);
    spdlog::trace("  energy = {}", energy);
    spdlog::trace("  emax = {}", emax);
    spdlog::trace("  s = {}", (point[0] + point[1]).M2());
    spdlog::trace("  s_min = {}", Smin());
#endif
}

double IntfSpectralMapper::GenerateWeight(const Event& event,
                                          std::vector<double> &rans) {
	const std::vector<FourVector>& point=event.Momentum();
    double pmin = point[0].E() - sqrt(pow(point[0].E(), 2) + 2 * point[0].E() * Constant::mN +
                                      Constant::mN2 - Smin());
    double pmax = point[0].E() + sqrt(pow(point[0].E(), 2) + 2 * point[0].E() * Constant::mN +
                                      Constant::mN2 - Smin());
    pmin = pmin < 0 ? 0 : pmin;
    pmax = pmax > 800 ? 800 : pmax;
    double dp = pmax - pmin;
    rans[0] = (point[HadronIdx()].P() - pmin) / dp;
    double cosT_max = (2 * point[0].E() * Constant::mN + Constant::mN2 - point[1].P2() - Smin()) /
                      (2 * point[0].E() * point[1].P());
    cosT_max = cosT_max > 1 ? 1 : cosT_max;
    double dCos = (cosT_max + 1);
    rans[1] = (point[HadronIdx()].CosTheta() + 1) / dCos;
    rans[2] = point[HadronIdx()].Phi() / dPhi;

    const double det =
        pow(point[0].E(), 2) + point[1].P2() + 2 * point[1].Vec3() * point[0].Vec3() + Smin();
    double emax = Constant::mN + point[0].E() - sqrt(det);
    emax = std::min(emax, Constant::mN - point[1].P());
    emax = emax > 400 ? 400 : emax;
    const double energy = Constant::mN - point[HadronIdx()].E();
    // if(energy < 0) return std::numeric_limits<double>::infinity();
    const double dE = emax;
    rans[3] = (energy + 1e-8) / emax;

    double spectator_momentum = point.back().P();
    rans[4] = spectator_momentum / dp2;
    rans[5] = (point.back().CosTheta() + 1) / dCos2;
    rans[6] = point.back().Phi() / dPhi;

    // rans[3] = (Constant::mN - point[HadronIdx()].E())/dE;

    // double cosT_max = (point[HadronIdx()].M2() +
    // 2*point[1].E()*point[HadronIdx()].E() -
    // Smin())/(2*point[HadronIdx()].P()*point[1].P()); cosT_max = cosT_max > 1
    // ? 1 : cosT_max; const double dCos = (cosT_max + 1);
    // rans[1] = (point[HadronIdx()].CosTheta() + 1) / dCos;

    double wgt =
        1.0 / point[1].P2() / dp / dCos / dPhi / dE / point.back().P2() / dp2 / dCos2 / dPhi;
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

DISSingleNucleonMapper::DISSingleNucleonMapper(const ProcessInfo& info,size_t idx): HadronicBeamMapper(info,idx) {
	SetMasses(info.Masses());
	// Hadron is always stationary in this model
	pHadron={ParticleInfo(info.m_hadronic.first[0]).Mass(),0,0,0};
}

void DISSingleNucleonMapper::GeneratePoint(Event& event,const std::vector<double>& rans) {
	std::vector<FourVector>& point=event.Momentum();
	// Generate initial quark state
	double x=xMin*pow(xMax/xMin,rans[0]);
	point[HadronIdx()]=getQuarkMomentum(point[0],pHadron,x);
}

double DISSingleNucleonMapper::GenerateWeight(const Event& event,std::vector<double>& rans) {
	const std::vector<FourVector>& point=event.Momentum();
	double x=getX(point[0],pHadron,point[HadronIdx()]);
	rans[0]=log(x/xMin)/log(xMax/xMin);
	return 1.0/(point[HadronIdx()].P2()*(x*log(xMax/xMin)));
}


DISNucleusMapper::DISNucleusMapper(const ProcessInfo& info,size_t idx): HadronicBeamMapper(info,idx) {
	SetMasses(info.Masses());
}

void DISNucleusMapper::GeneratePoint(Event& event,const std::vector<double>& rans) {
	std::vector<FourVector>& point=event.Momentum();
	// Generate inital nucleon state
	double radical =
		pow(point[0].E(), 2) + 2 * point[0].E() * Constant::mN + Constant::mN2 - Smin();
	if(radical < 0) radical = 0;
	double pmin = point[0].E() - sqrt(radical);
	double pmax = point[0].E() + sqrt(radical);
	pmin = pmin < 0 ? 0 : pmin;
	pmax = pmax > 800 ? 800 : pmax;
	double dp = pmax - pmin;
	const double mom = dp * rans[0] + pmin;
	double cosT_max = (2 * point[0].E() * Constant::mN + Constant::mN2 - mom * mom - Smin()) /
					  (2 * point[0].E() * mom);
	cosT_max = cosT_max > 1 ? 1 : cosT_max < -1 ? -1 : cosT_max;
	const double cosT = (cosT_max + 1) * rans[1] - 1;
	const double sinT = sqrt(1 - cosT * cosT);
	const double phi = dPhi * rans[2];
	ThreeVector pmom = {mom * sinT * cos(phi), mom * sinT * sin(phi), mom * cosT};

	const double det = pow(point[0].E(), 2) + mom * mom + 2 * pmom * point[0].Vec3() + Smin();
	double emax = Constant::mN + point[0].E() - sqrt(det);
	emax = std::min(emax, Constant::mN - mom);
	emax = emax > 400 ? 400 : emax;
	const double energy = emax * rans[3] - 1e-8;

	// pHadron
	point[HadronIdx()]={Constant::mN - energy, mom * sinT * cos(phi), mom * sinT * sin(phi), mom * cosT};

	double x=xMin*pow(xMax/xMin,rans[4]);

	// pQuark
	point[HadronIdx()+1]=getQuarkMomentum(point[0],point[HadronIdx()],x);
}

double DISNucleusMapper::GenerateWeight(const Event& event,std::vector<double>& rans) {
	const std::vector<FourVector>& point=event.Momentum();
	double x=getX(point[0],point[HadronIdx()],point[HadronIdx()+1]);

	double radical =
	    pow(point[0].E(), 2) + 2 * point[0].E() * Constant::mN + Constant::mN2 - Smin();
	if(radical < 0) radical = 0;
	double pmin = point[0].E() - sqrt(radical);
	double pmax = point[0].E() + sqrt(radical);
	pmin = pmin < 0 ? 0 : pmin;
	pmax = pmax > 800 ? 800 : pmax;
	double dp = pmax - pmin;
	rans[0] = (point[HadronIdx()].P() - pmin) / dp;
	double cosT_max = (2 * point[0].E() * Constant::mN + Constant::mN2 - point[HadronIdx()].P2() - Smin()) /
	                  (2 * point[0].E() * point[HadronIdx()].P());
	cosT_max = cosT_max > 1 ? 1 : cosT_max < -1 ? -1 : cosT_max;
	double dCos = (cosT_max + 1);
	rans[1] = (point[HadronIdx()].CosTheta() + 1) / dCos;
	rans[2] = point[HadronIdx()].Phi() / dPhi;

	const double det =
	    pow(point[0].E(), 2) + point[HadronIdx()].P2() + 2 * point[HadronIdx()].Vec3() * point[0].Vec3() + Smin();
	double emax = Constant::mN + point[0].E() - sqrt(det);
	emax = std::min(emax, Constant::mN - point[HadronIdx()].P());
	emax = emax > 400 ? 400 : emax;
	const double energy = Constant::mN - point[HadronIdx()].E();
	const double dE = emax;
	rans[3] = (energy + 1e-8) / emax;

	rans[4]=log(x/xMin)/log(xMax/xMin);

    return 1.0/(point[HadronIdx()+1].P2()*dp*dCos*dPhi*dE*(x*log(xMax/xMin)));
}