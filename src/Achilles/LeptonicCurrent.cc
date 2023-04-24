#include "Achilles/LeptonicCurrent.hh"
#include "Achilles/Constants.hh"
#include "Achilles/FormFactor.hh"
#include "Achilles/FourVector.hh"
#include "Achilles/ProcessInfo.hh"
#include "Achilles/Spinor.hh"

using achilles::LeptonicCurrent;

void LeptonicCurrent::Initialize(const ProcessInfo &process) {
    using namespace achilles::Constant;
    const std::complex<double> i(0, 1);
    // Determine process
    bool init_neutrino = ParticleInfo(process.m_leptonic.first).IsNeutrino();
    bool neutral_current = NeutralCurrent(process.m_leptonic.first, process.m_leptonic.second[0]);
    bool charged_current =
        ChargedCurrent(init_neutrino, process.m_leptonic.first, process.m_leptonic.second[0]);
    if(!neutral_current && !charged_current)
        throw std::runtime_error("HardScattering: Invalid process");

    // TODO: Define couplings correctly
    if(charged_current) {
        pid = init_neutrino ? (process.m_leptonic.first.AsInt() < 0 ? 24 : -24)
                            : (process.m_leptonic.first.AsInt() < 0 ? -24 : 24);
        coupl_right = 0;
        coupl_left = ee * i / (sw * sqrt(2));
        mass = Constant::MW;
        width = Constant::GAMW;
    } else if(neutral_current) {
        if(init_neutrino) {
            coupl_left = (cw * ee * i) / (2 * sw) + (ee * i * sw) / (2 * cw);
            coupl_right = 0;
            pid = 23;
            mass = Constant::MZ;
            width = Constant::GAMZ;
        } else {
            coupl_right = -ee * i;
            coupl_left = coupl_right;
            pid = 22;
        }
    }
    anti = process.m_leptonic.first.AsInt() < 0;
}

bool LeptonicCurrent::NeutralCurrent(achilles::PID initial, achilles::PID final) const {
    return initial == final;
}

bool LeptonicCurrent::ChargedCurrent(bool neutrino, achilles::PID initial,
                                     achilles::PID final) const {
    return initial.AsInt() - (2 * neutrino - 1) == final.AsInt();
}

achilles::FFDictionary LeptonicCurrent::GetFormFactor() {
    FFDictionary results;
    static constexpr std::complex<double> i(0, 1);
    using namespace achilles::Constant;
    // TODO: Double check form factors
    if(pid == 24) {
        const std::complex<double> coupl = ee * i / (sw * sqrt(2) * 2);
        results[{PID::proton(), pid}] = {{FormFactorInfo::Type::F1p, coupl},
                                         {FormFactorInfo::Type::F1n, -coupl},
                                         {FormFactorInfo::Type::F2p, coupl},
                                         {FormFactorInfo::Type::F2n, -coupl},
                                         {FormFactorInfo::Type::FA, coupl}};
        results[{PID::neutron(), pid}] = {};
        results[{PID::carbon(), pid}] = {};
    } else if(pid == -24) {
        const std::complex<double> coupl = ee * i / (sw * sqrt(2) * 2);
        results[{PID::neutron(), pid}] = {{FormFactorInfo::Type::F1p, coupl},
                                          {FormFactorInfo::Type::F1n, -coupl},
                                          {FormFactorInfo::Type::F2p, coupl},
                                          {FormFactorInfo::Type::F2n, -coupl},
                                          {FormFactorInfo::Type::FA, coupl}};
        results[{PID::proton(), pid}] = {};
        results[{PID::carbon(), pid}] = {};
    } else if(pid == 23) {
        const std::complex<double> coupl1 = cw * ee * i / (2 * sw) - ee * i * sw / (2 * cw);
        const std::complex<double> coupl2 = -(cw * ee * i / (2 * sw));
        results[{PID::proton(), pid}] = {{FormFactorInfo::Type::F1p, coupl1},
                                         {FormFactorInfo::Type::F1n, coupl2},
                                         {FormFactorInfo::Type::F2p, coupl1},
                                         {FormFactorInfo::Type::F2n, coupl2},
                                         {FormFactorInfo::Type::FA, coupl2}};
        results[{PID::neutron(), pid}] = {{FormFactorInfo::Type::F1n, coupl1},
                                          {FormFactorInfo::Type::F1p, coupl2},
                                          {FormFactorInfo::Type::F2n, coupl1},
                                          {FormFactorInfo::Type::F2p, coupl2},
                                          {FormFactorInfo::Type::FA, coupl2}};
        results[{PID::carbon(), pid}] = {};
    } else if(pid == 22) {
        const std::complex<double> coupl = i * ee;
        results[{PID::proton(), pid}] = {{FormFactorInfo::Type::F1p, coupl},
                                         {FormFactorInfo::Type::F2p, coupl}};
        results[{PID::neutron(), pid}] = {{FormFactorInfo::Type::F1n, coupl},
                                          {FormFactorInfo::Type::F2n, coupl}};
        results[{PID::carbon(), pid}] = {{FormFactorInfo::Type::FCoh, 6.0 * coupl}};
    } else {
        throw std::runtime_error("LeptonicCurrent: Invalid probe");
    }

    return results;
}

achilles::Currents LeptonicCurrent::CalcCurrents(const FourVector &p_in,
                                                 const FourVector &p_out) const {
    Currents currents;

    // Setup spinors
    FourVector pU, pUBar;
    if(anti) {
        pUBar = -p_in;
        pU = p_out;
    } else {
        pU = -p_in;
        pUBar = p_out;
    }
    std::array<Spinor, 2> ubar, u;
    ubar[0] = UBarSpinor(-1, pUBar);
    ubar[1] = UBarSpinor(1, pUBar);
    u[0] = USpinor(-1, pU);
    u[1] = USpinor(1, pU);

    // Calculate currents
    Current result;
    double q2 = (p_in - p_out).M2();
    std::complex<double> prop =
        std::complex<double>(0, 1) / (q2 - mass * mass - std::complex<double>(0, 1) * mass * width);
    spdlog::trace("Calculating Current for {}", pid);
    for(size_t i = 0; i < 2; ++i) {
        for(size_t j = 0; j < 2; ++j) {
            VCurrent subcur;
            for(size_t mu = 0; mu < 4; ++mu) {
                subcur[mu] = ubar[i] *
                             (coupl_left * SpinMatrix::GammaMu(mu) * SpinMatrix::PL() +
                              coupl_right * SpinMatrix::GammaMu(mu) * SpinMatrix::PR()) *
                             u[j] * prop;
                spdlog::trace("Current[{}][{}] = {}", 2 * i + j, mu, subcur[mu]);
            }
            result.push_back(subcur);
        }
    }
    currents[pid] = result;

    return currents;
}
