#include "nuchic/Cuts.hh"
#include "nuchic/FourVector.hh"

bool nuchic::Cut::operator()(const FourVector &vec) const {
    for(const auto &cut : m_cuts) {
        switch(cut.first) {
            case CutType::Q2:
                if(!Q2Cut(vec, cut.second)) return false;
                break;
            case CutType::Energy:
                if(!EnergyCut(vec, cut.second)) return false;
                break;
            case CutType::Momentum:
                if(!MomentumCut(vec, cut.second)) return false;
                break;
            case CutType::InvariantMass:
                if(!InvariantMassCut(vec, cut.second)) return false;
                break;
            case CutType::TransverseMomentum:
                if(!TransverseMomentumCut(vec, cut.second)) return false;
                break;
            case CutType::AngleTheta:
                if(!AngleThetaCut(vec, cut.second)) return false;
                break;
            case CutType::AnglePhi:
                if(!AnglePhiCut(vec, cut.second)) return false;
                break;
            case CutType::ETheta2:
                if(!ETheta2Cut(vec, cut.second)) return false;
                break;
        }
    }
    return true;
}

bool nuchic::Cut::MakeCut(double val, const cut_range &cutRange) const {
    bool result = false;

    for(const auto &cut : cutRange) {
        result |= (val > cut.first && val < cut.second);
    }

    return result;
}

bool nuchic::Cut::Q2Cut(const FourVector &/*vec*/, const cut_range &/*cutRange*/) const {
    // TODO: Implement Q2 cut
    return true;
}

bool nuchic::Cut::EnergyCut(const FourVector &vec, const cut_range &cutRange) const {
    return MakeCut(vec.E(), cutRange);
}

bool nuchic::Cut::MomentumCut(const FourVector &vec, const cut_range &cutRange) const {
    return MakeCut(vec.P(), cutRange);
}

bool nuchic::Cut::InvariantMassCut(const FourVector &vec, const cut_range &cutRange) const {
    return MakeCut(vec.M(), cutRange);
}

bool nuchic::Cut::TransverseMomentumCut(const FourVector &vec, const cut_range &cutRange) const {
    return MakeCut(vec.Pt(), cutRange);
}

bool nuchic::Cut::AngleThetaCut(const FourVector &vec, const cut_range &cutRange) const {
    return MakeCut(vec.Theta(), cutRange);
}

bool nuchic::Cut::AnglePhiCut(const FourVector &vec, const cut_range &cutRange) const {
    return MakeCut(vec.Phi(), cutRange);
}

bool nuchic::Cut::ETheta2Cut(const FourVector &vec, const cut_range &cutRange) const {
    return MakeCut(vec.E()*pow(vec.Theta(),2), cutRange);
}

void nuchic::Cut::Add(const CutType &type, double cutVal) {
    if(m_cuts.find(type) != m_cuts.end()) {
        auto errorMsg = fmt::format("Cut: May only have one cut of each type except Custom. Found additional cut of type: {}", type);
        throw std::runtime_error(errorMsg);
    }

    m_cuts[type] = {{cutVal, std::numeric_limits<double>::infinity()}};
}

void nuchic::Cut::Add(const CutType &type, cut_range cutVal) {
    if(m_cuts.find(type) != m_cuts.end()) {
        auto errorMsg = fmt::format("Cut: May only have one cut of each type except Custom. Found additional cut of type: {}", type);
        throw std::runtime_error(errorMsg);
    }

    auto& lastCut = cutVal[cutVal.size() - 1];
    if(lastCut.first > 0 && lastCut.second < 0) {
        lastCut.second = std::numeric_limits<double>::infinity(); 
    }
    m_cuts[type] = cutVal;
}

void nuchic::Cut::Add(const std::function<double(const nuchic::FourVector&)> &/*func*/, double /*cutVal*/) {
    // TODO: Implement handling of custom cuts
}

void nuchic::Cut::Add(const std::function<double(const nuchic::FourVector&)> &/*func*/, cut_range /*cutVal*/) {
    // TODO: Implement handling of custom cuts
}
