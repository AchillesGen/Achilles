#ifdef ACHILLES_EVENT_DETAILS
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_TRACE
#endif

#include "Plugins/Sherpa/Channels3.hh"
#include "ATOOLS/Phys/Flavour.H"
#include "Achilles/Utilities.hh"
#include "COMIX/Main/Single_Process.H"
#include "Plugins/Sherpa/PrintVec.hh"
#include <cmath>

using namespace PHASIC;
using namespace ATOOLS;

double GenChannel::SCut(size_t id) {
    if(id & 3) { id = (1 << m_n) - 1 - id; }
    double result = 0;
    for(size_t i = 0; i < m_n; ++i) {
        if(achilles::SetBit(id, i)) { result += sqrt(m_s[i]); }
    }
    return result * result;
}

size_t GenChannel::SId(size_t id) {
    return (id & 3) == 3 ? (1 << m_n) - 1 - id : id;
}

void GenChannel::FillMomenta(ChannelNode *node) {
    size_t aid = (1 << m_n) - 1 - node->m_idx;
    m_p[node->m_idx] = Vec4D{};
    for(size_t i = 0; i < m_n; ++i) {
        if(achilles::SetBit(node->m_idx, i)) { m_p[node->m_idx] += m_p[(1 << i)]; }
    }
    if(!achilles::IsPower2(node->m_left->m_idx)) FillMomenta(node->m_left.get());
    if(!achilles::IsPower2(node->m_right->m_idx)) FillMomenta(node->m_right.get());
}

std::string GenChannel::PrintPoint(ChannelNode *cur, size_t lid, size_t depth) const {
    std::string result = "";
    if(depth == m_n - 2) return result;
    size_t aid = cur->m_idx;
    size_t bid = cur->m_left->m_idx;
    size_t cid = cur->m_right->m_idx;
    if(bid == lid)
        std::swap(aid, bid);
    else if(cid == lid)
        std::swap(aid, cid);
    if((cid & (lid | m_rid)) == (lid | m_rid) || (aid & m_rid && bid & m_rid)) {
        std::swap(bid, cid);
        std::swap(cur->m_left, cur->m_right);
    }
    if(cid == m_rid) {
        if(!achilles::IsPower2(bid)) {
            result += PrintSPoint(cur->m_left.get(), bid);
            return result;
        } else {
            return result;
        }
    }
    result += fmt::format("TChannel({} ({}), {} ({}), {} ({})),", aid, cur->m_pid, bid,
                          cur->m_left->m_pid, cid, cur->m_right->m_pid);
    result += PrintPoint(cur->m_right.get(), cid, depth + 1);
    return result;
}

std::string GenChannel::PrintSPoint(ChannelNode *node, size_t id) const {
    size_t cid = node->m_idx;
    size_t aid = node->m_left->m_idx;
    size_t bid = node->m_right->m_idx;
    std::string result = fmt::format("SChannel({} ({}), {} ({}), {} ({})),", cid, node->m_pid, aid,
                                     node->m_left->m_pid, bid, node->m_right->m_pid);
    if(!achilles::IsPower2(aid)) result += PrintSPoint(node->m_left.get(), aid);
    if(!achilles::IsPower2(bid)) result += PrintSPoint(node->m_right.get(), bid);
    return result;
}

void GenChannel::GeneratePoint(std::vector<ATOOLS::Vec4D> &p, const std::vector<double> &rans) {
    iran = 0;
    m_p[1] = p[0];
    m_p[2] = p[1];
    auto lid = static_cast<size_t>((1 << m_n) - 2);
    m_p[lid] = p[0];
    BuildPoint(m_nodes.get(), lid, rans, 0);
    for(size_t i = 2; i < m_n; ++i) { p[i] = m_p[1 << i]; }
#ifdef ACHILLES_EVENT_DETAILS
    Mapper::Print(__PRETTY_FUNCTION__, p, rans);
#endif
}

void GenChannel::BuildPoint(ChannelNode *cur, size_t lid, const std::vector<double> &rans,
                            size_t depth) {
    if(depth == m_n - 2) { return; }
    size_t aid = cur->m_idx;
    size_t bid = cur->m_left->m_idx;
    size_t cid = cur->m_right->m_idx;
    if(bid == lid)
        std::swap(aid, bid);
    else if(cid == lid)
        std::swap(aid, cid);
    if((cid & (lid | m_rid)) == (lid | m_rid) || (aid & m_rid && bid & m_rid)) {
        std::swap(bid, cid);
        std::swap(cur->m_left, cur->m_right);
    }
    if(cid == m_rid) {
        if(!achilles::IsPower2(bid)) {
            BuildSPoint(cur->m_left.get(), bid, rans);
            return;
        } else {
            return;
        }
    }
    TChannelMomenta(cur, aid, bid, cid, rans);
    BuildPoint(cur->m_right.get(), cid, rans, depth + 1);
}

void GenChannel::BuildSPoint(ChannelNode *node, size_t id, const std::vector<double> &rans) {
    size_t cid = node->m_idx;
    size_t aid = node->m_left->m_idx;
    size_t bid = node->m_right->m_idx;
    SChannelMomenta(node, aid, bid, cid, rans);
    if(!achilles::IsPower2(aid)) BuildSPoint(node->m_left.get(), aid, rans);
    if(!achilles::IsPower2(bid)) BuildSPoint(node->m_right.get(), bid, rans);
}

double GenChannel::BuildWeight(ChannelNode *cur, size_t lid, std::vector<double> &rans,
                               size_t depth) {
    double wgt = 1.0;
    if(depth == m_n - 2) { return wgt; }
    size_t aid = cur->m_idx;
    size_t bid = cur->m_left->m_idx;
    size_t cid = cur->m_right->m_idx;
    if(bid == lid)
        std::swap(aid, bid);
    else if(cid == lid)
        std::swap(aid, cid);
    if((cid & (lid | m_rid)) == (lid | m_rid) || (aid & m_rid && bid & m_rid)) {
        std::swap(bid, cid);
        std::swap(cur->m_left, cur->m_right);
    }
    if(cid == m_rid) {
        if(!achilles::IsPower2(bid)) {
            wgt *= BuildSWeight(cur->m_left.get(), bid, rans);
            return wgt;
        } else {
            return wgt;
        }
    }
    wgt *= TChannelWeight(cur, aid, bid, cid, rans);
    wgt *= BuildWeight(cur->m_right.get(), cid, rans, depth + 1);
    return wgt;
}

double GenChannel::BuildSWeight(ChannelNode *node, size_t id, std::vector<double> &rans) {
    double wgt = 1.0;
    size_t cid = node->m_idx;
    size_t aid = node->m_left->m_idx;
    size_t bid = node->m_right->m_idx;
    wgt *= SChannelWeight(node, aid, bid, cid, rans);
    if(!achilles::IsPower2(aid)) wgt *= BuildSWeight(node->m_left.get(), aid, rans);
    if(!achilles::IsPower2(bid)) wgt *= BuildSWeight(node->m_right.get(), bid, rans);
    return wgt;
}

double GenChannel::PropMomenta(ChannelNode *node, size_t id, double smin, double smax, double ran) {
    auto foundNode = LocateNode(node, id);
    if(!foundNode) { return CE.ThresholdMomenta(m_salpha, 0.01, smin, smax, ran); }
    auto flav = Flavour((kf_code)(foundNode->m_pid));
    if(flav.Mass() == 0) {
        return CE.MasslessPropMomenta(m_salpha, smin, smax, ran);
    } else {
        return CE.MassivePropMomenta(flav.Mass(), flav.Width(), 1, smin, smax, ran);
    }
}

double GenChannel::PropWeight(ChannelNode *node, size_t id, double smin, double smax, double s,
                              double &ran) {
    auto foundNode = LocateNode(node, id);
    if(!foundNode) {
        double wgt = CE.ThresholdWeight(m_salpha, 0.01, smin, smax, s, ran);
        return wgt;
    }
    auto flav = Flavour((kf_code)(foundNode->m_pid));
    double wgt = 1.0;
    if(flav.Mass() == 0) {
        wgt = CE.MasslessPropWeight(m_salpha, smin, smax, s, ran);
    } else {
        wgt = CE.MassivePropWeight(flav.Mass(), flav.Width(), 1, smin, smax, s, ran);
    }

    return wgt;
}

void GenChannel::TChannelMomenta(ChannelNode *node, size_t aid, size_t bid, size_t cid,
                                 const std::vector<double> &rans) {
    size_t pid = aid - (m_rid + bid);
    const std::string name = "TChannelMomenta(" + std::to_string(aid) + ", " + std::to_string(bid) +
                             ", " + std::to_string(cid) + ", " + std::to_string(pid) + ")";
    double se = SCut(bid), sp = SCut(pid);
    double rtsmax = (m_p[aid] + m_p[m_rid]).Mass();
    if(!achilles::IsPower2(bid)) {
        double smin = se, smax = sqr(rtsmax - sqrt(sp));
        se = PropMomenta(node, bid, smin, smax, rans[iran++]);
    }
    if(!achilles::IsPower2(pid)) {
        double smin = sp, smax = sqr(rtsmax - sqrt(se));
        sp = PropMomenta(node, pid, smin, smax, rans[iran++]);
    }
    double mass = ATOOLS::Flavour((kf_code)(LocateNode(node, pid + m_rid)->m_pid)).Mass();
    CE.TChannelMomenta(m_p[aid], m_p[m_rid], m_p[bid], m_p[pid], se, sp, mass, m_alpha, m_ctmax,
                       m_ctmin, m_amct, 0, rans[iran], rans[iran + 1]);
    iran += 2;
    m_p[cid] = m_p[aid] - m_p[bid];
    SPDLOG_TRACE("{}", name);
    SPDLOG_TRACE("  m_p[{}] = {}, m = {}", aid, m_p[aid], m_p[aid].Mass());
    SPDLOG_TRACE("  m_p[{}] = {}, m = {}", m_rid, m_p[m_rid], m_p[m_rid].Mass());
    SPDLOG_TRACE("  m_p[{}] = {}, m = {}", bid, m_p[bid], m_p[bid].Mass());
    SPDLOG_TRACE("  m_p[{}] = {}, m = {}", pid, m_p[pid], m_p[pid].Mass());
    SPDLOG_TRACE("  se = {}, sp = {}", se, sp);
    SPDLOG_TRACE("  m[0]^2 = {}, m[1]^2 = {}, m[2]^2 = {}, m[3]^2 = {}", m_s[0], m_s[1], m_s[2],
                 m_s[3]);
    SPDLOG_TRACE("  iran = {}", iran);
}

double GenChannel::TChannelWeight(ChannelNode *node, size_t aid, size_t bid, size_t cid,
                                  std::vector<double> &rans) {
    size_t pid = aid - (m_rid + bid);
    double wgt = 1.0;
    // aid = (1 << m_n) - 1 - aid;
    m_p[pid] = m_p[aid] + m_p[m_rid] - m_p[bid];
    double se = SCut(bid), sp = SCut(pid);
    double rtsmax = (m_p[aid] + m_p[m_rid]).Mass();
    const std::string name = "TChannelWeight(" + std::to_string(aid) + ", " + std::to_string(bid) +
                             ", " + std::to_string(cid) + ", " + std::to_string(pid) + ")";
    SPDLOG_TRACE("{}", name);
    if(!achilles::IsPower2(bid)) {
        double smin = se, smax = sqr(rtsmax - sqrt(sp));
        wgt *= PropWeight(node, bid, smin, smax, se = m_p[bid].Abs2(), rans[iran++]);
        SPDLOG_TRACE("  smin = {}, smax = {}", smin, smax);
    }
    if(!achilles::IsPower2(pid)) {
        double smin = sp, smax = sqr(rtsmax - sqrt(se));
        wgt *= PropWeight(node, pid, smin, smax, sp = m_p[pid].Abs2(), rans[iran++]);
        SPDLOG_TRACE("  smin = {}, smax = {}", smin, smax);
    }
    double mass = ATOOLS::Flavour((kf_code)(LocateNode(node, pid + m_rid)->m_pid)).Mass();
    wgt *= CE.TChannelWeight(m_p[aid], m_p[m_rid], m_p[bid], m_p[pid], mass, m_alpha, m_ctmax,
                             m_ctmin, m_amct, 0, rans[iran], rans[iran + 1]);
    iran += 2;
    m_p[cid] = m_p[aid] - m_p[bid];
    SPDLOG_TRACE("  m_p[{}] = {}", aid, m_p[aid]);
    SPDLOG_TRACE("  m_p[{}] = {}", m_rid, m_p[m_rid]);
    SPDLOG_TRACE("  m_p[{}] = {}", bid, m_p[bid]);
    SPDLOG_TRACE("  m_p[{}] = {}", pid, m_p[pid]);
    SPDLOG_TRACE("  mass = {}", mass);
    SPDLOG_TRACE("  se = {}, sp = {}", se, sp);
    SPDLOG_TRACE("  iran = {}", iran);
    SPDLOG_TRACE("  wgt = {}", wgt);
    return wgt;
}

void GenChannel::SChannelMomenta(ChannelNode *node, size_t aid, size_t bid, size_t cid,
                                 const std::vector<double> &rans) {
    size_t lid = SId(aid), rid = SId(bid);
    const std::string name = "SChannelMomenta(" + std::to_string(cid) + ", " + std::to_string(aid) +
                             ", " + std::to_string(bid) + ")";
    double rts = m_p[cid].Mass(), sl = SCut(lid), sr = SCut(rid);
    if(!achilles::IsPower2(lid)) {
        double smin = sl, smax = sqr(rts - sqrt(sr));
        SPDLOG_TRACE("smin = {}, smax = {}, rts = {}, sr = {}", smin, smax, rts, sr);
        sl = PropMomenta(node, lid, smin, smax, rans[iran++]);
    }
    if(!achilles::IsPower2(rid)) {
        double smin = sr, smax = sqr(rts - sqrt(sl));
        SPDLOG_TRACE("smin = {}, smax = {}, rts = {}, sl = {}", smin, smax, rts, sl);
        sr = PropMomenta(node, rid, smin, smax, rans[iran++]);
    }
    CE.Isotropic2Momenta(m_p[cid], sl, sr, m_p[lid], m_p[rid], rans[iran], rans[iran + 1], m_ctmin,
                         m_ctmax);
    iran += 2;
    m_p[(1 << m_n) - 1 - aid] = m_p[aid];
    m_p[(1 << m_n) - 1 - bid] = m_p[bid];
    SPDLOG_TRACE("{}", name);
    SPDLOG_TRACE("  m_p[{}] ({}) = {}, m = {}", cid, node->m_pid, m_p[cid], m_p[cid].Mass());
    SPDLOG_TRACE("  m_p[{}] ({}) = {}, m = {}", lid, node->m_left->m_pid, m_p[lid],
                 m_p[lid].Mass());
    SPDLOG_TRACE("  m_p[{}] ({}) = {}, m = {}", rid, node->m_right->m_pid, m_p[rid],
                 m_p[rid].Mass());
    SPDLOG_TRACE("  sl = {}, sr = {}", sl, sr);
    SPDLOG_TRACE("  iran = {}", iran);
}

double GenChannel::SChannelWeight(ChannelNode *node, size_t aid, size_t bid, size_t cid,
                                  std::vector<double> &rans) {
    double wgt = 1.0;
    size_t lid = SId(aid), rid = SId(bid);
    double rts = m_p[cid].Mass(), sl = SCut(lid), sr = SCut(rid);
    if(!achilles::IsPower2(lid)) {
        double smin = sl, smax = sqr(rts - sqrt(sr));
        wgt *= PropWeight(node, lid, smin, smax, sl = m_p[lid].Abs2(), rans[iran++]);
    }
    if(!achilles::IsPower2(rid)) {
        double smin = sr, smax = sqr(rts - sqrt(sl));
        wgt *= PropWeight(node, rid, smin, smax, sr = m_p[rid].Abs2(), rans[iran++]);
    }
    wgt *= CE.Isotropic2Weight(m_p[lid], m_p[rid], rans[iran], rans[iran + 1], m_ctmin, m_ctmax);
    iran += 2;
    const std::string name = "SChannelWeight(" + std::to_string(cid) + ", " + std::to_string(aid) +
                             ", " + std::to_string(bid) + ")";
    SPDLOG_TRACE("{}", name);
    SPDLOG_TRACE("  m_p[{}] ({}) = {}", cid, node->m_pid, m_p[cid]);
    SPDLOG_TRACE("  m_p[{}] ({}) = {}", lid, node->m_left->m_pid, m_p[lid]);
    SPDLOG_TRACE("  m_p[{}] ({}) = {}", rid, node->m_right->m_pid, m_p[rid]);
    SPDLOG_TRACE("  sl = {}, sr = {}", sl, sr);
    SPDLOG_TRACE("  iran = {}", iran);
    SPDLOG_TRACE("  wgt = {}", wgt);
    return wgt;
}

ChannelNode *GenChannel::LocateNode(ChannelNode *node, size_t id) {
    if(node->m_idx == id) return node;
    if(node->m_left) {
        auto result = LocateNode(node->m_left.get(), id);
        if(!result && node->m_right) { result = LocateNode(node->m_right.get(), id); }
        return result;
    }
    return nullptr;
}

bool GenChannel::InitializeChannel(std::shared_ptr<ChannelNode> node) {
    m_nodes = node;
    return true;
}

void Channels3::GenerateNuclearPoint(std::vector<Vec4D> &p, const std::vector<double> &ran) const {
    double s234_max = sqr((p[0] + p[1]).Mass() - sqrt(s5));
    double s234_min = pow(sqrt(s3) + sqrt(s4) + sqrt(s2), 2);
    double s234 = CE.MasslessPropMomenta(0.5, s234_min, s234_max, ran[5]);
    Vec4D p234;
    CE.TChannelMomenta(p[0], p[1], p[5], p234, s5, s234, 0., m_alpha, m_ctmax, m_ctmin, m_amct, 0,
                       ran[6], ran[7]);
}

double Channels3::GenerateNuclearWeight(const std::vector<Vec4D> &p,
                                        std::vector<double> &rans) const {
    double wt = 1.;

    double s234_max = sqr((p[0] + p[1]).Mass() - sqrt(s5));
    double s234_min = pow(sqrt(s3) + sqrt(s4) + sqrt(s2), 2);
    wt *=
        CE.MasslessPropWeight(0.5, s234_min, s234_max, dabs((p[2] + p[3] + p[4]).Abs2()), rans[5]);
    wt *= CE.TChannelWeight(p[0], p[1], p[5], p[2] + p[3] + p[4], 0., m_alpha, m_ctmax, m_ctmin,
                            m_amct, 0, rans[6], rans[7]);

    return wt;
}

void C3_0::GeneratePoint(std::vector<Vec4D> &p, const std::vector<double> &ran) {
    GenerateNuclearPoint(p, ran);
    double s34_max = sqr((p[0] + p[1] - p[5]).Mass() - sqrt(s2));
    // TODO: set minimum from cuts
    // double s34_min = cuts->GetscutAmegic(std::string("34"));
    double s34_min = std::max(pow(sqrt(s3) + sqrt(s4), 2), 1E-8);
    Vec4D p34;
    double s34 = CE.MasslessPropMomenta(m_salpha, s34_min, s34_max, ran[0]);
    CE.Isotropic2Momenta(p[0] + p[1] - p[5], s2, s34, p[2], p34, ran[1], ran[2]);
    CE.Isotropic2Momenta(p34, s3, s4, p[3], p[4], ran[3], ran[4]);
    // Mapper<Vec4D>::Print(__PRETTY_FUNCTION__, p, ran);
}

double C3_0::GenerateWeight(const std::vector<Vec4D> &p, std::vector<double> &rans) {
    double wt = GenerateNuclearWeight(p, rans);
    double s34_max = sqr((p[0] + p[1] - p[5]).Mass() - sqrt(s2));
    // TODO: set minimum from cuts
    // double s34_min = cuts->GetscutAmegic(std::string("34"));
    double s34_min = std::max(pow(sqrt(s3) + sqrt(s4), 2), 1E-8);
    wt *= CE.MasslessPropWeight(m_salpha, s34_min, s34_max, dabs((p[3] + p[4]).Abs2()), rans[0]);
    wt *= CE.Isotropic2Weight(p[2], p[3] + p[4], rans[1], rans[2], -1, 1);
    wt *= CE.Isotropic2Weight(p[3], p[4], rans[3], rans[4], -1, 1);
    if(wt != 0.) wt = 1.0 / wt / pow(2. * M_PI, 3 * 4. - 4.);
    // Mapper<Vec4D>::Print(__PRETTY_FUNCTION__, p, rans);
    return 1 / wt;
}

void C3_1::GeneratePoint(std::vector<Vec4D> &p, const std::vector<double> &ran) {
    GenerateNuclearPoint(p, ran);
    double s34_max = sqr((p[0] + p[1] - p[5]).Mass() - sqrt(s2));
    // TODO: set minimum from cuts
    // double s34_min = cuts->GetscutAmegic(std::string("34"));
    double s34_min = std::max(pow(sqrt(s3) + sqrt(s4), 2), 1E-8);
    Flavour fl34 = Flavour((kf_code)(23));
    Vec4D p34;
    double s34 = CE.MassivePropMomenta(fl34.Mass(), fl34.Width(), 1, s34_min, s34_max, ran[0]);
    CE.Isotropic2Momenta(p[0] + p[1] - p[5], s2, s34, p[2], p34, ran[1], ran[2]);
    CE.Isotropic2Momenta(p34, s3, s4, p[3], p[4], ran[3], ran[4]);
    // Mapper<Vec4D>::Print(__PRETTY_FUNCTION__, p, ran);
}

double C3_1::GenerateWeight(const std::vector<Vec4D> &p, std::vector<double> &rans) {
    double wt = GenerateNuclearWeight(p, rans);
    double s34_max = sqr((p[0] + p[1] - p[5]).Mass() - sqrt(s2));
    // TODO: set minimum from cuts
    // double s34_min = cuts->GetscutAmegic(std::string("34"));
    double s34_min = std::max(pow(sqrt(s3) + sqrt(s4), 2), 1E-8);
    Flavour fl34 = Flavour((kf_code)(23));
    wt *= CE.MassivePropWeight(fl34.Mass(), fl34.Width(), 1, s34_min, s34_max,
                               dabs((p[3] + p[4]).Abs2()), rans[0]);
    wt *= CE.Isotropic2Weight(p[2], p[3] + p[4], rans[1], rans[2], -1, 1);
    wt *= CE.Isotropic2Weight(p[3], p[4], rans[3], rans[4], -1, 1);
    if(wt != 0.) wt = 1.0 / wt / pow(2. * M_PI, 3 * 4. - 4.);
    // Mapper<Vec4D>::Print(__PRETTY_FUNCTION__, p, rans);
    return 1 / wt;
}

void C3_2::GeneratePoint(std::vector<Vec4D> &p, const std::vector<double> &ran) {
    GenerateNuclearPoint(p, ran);
    double s34_max = sqr((p[0] + p[1] - p[5]).Mass() - sqrt(s2));
    // TODO: set minimum from cuts
    // double s34_min = cuts->GetscutAmegic(std::string("34"));
    double s34_min = std::max(pow(sqrt(s3) + sqrt(s4), 2), 1E-8);
    Vec4D p34;
    double s34 = CE.MasslessPropMomenta(m_salpha, s34_min, s34_max, ran[0]);
    CE.TChannelMomenta(p[0] - p[5], p[1], p[2], p34, s2, s34, 0., m_alpha, m_ctmax, m_ctmin, m_amct,
                       0, ran[1], ran[2]);
    CE.TChannelMomenta(p[0] - p[5], p[1] - p[2], p[4], p[3], s4, s3, 0., m_alpha, 1., -1., m_amct,
                       0, ran[3], ran[4]);
    // Mapper<Vec4D>::Print(__PRETTY_FUNCTION__, p, ran);
}

double C3_2::GenerateWeight(const std::vector<Vec4D> &p, std::vector<double> &rans) {
    double wt = GenerateNuclearWeight(p, rans);
    double s34_max = sqr((p[0] + p[1] - p[5]).Mass() - sqrt(s2));
    // TODO: set minimum from cuts
    // double s34_min = cuts->GetscutAmegic(std::string("34"));
    double s34_min = std::max(pow(sqrt(s3) + sqrt(s4), 2), 1E-8);
    wt *= CE.MasslessPropWeight(m_salpha, s34_min, s34_max, dabs((p[3] + p[4]).Abs2()), rans[0]);
    wt *= CE.TChannelWeight(p[0] - p[5], p[1], p[2], p[3] + p[4], 0., m_alpha, m_ctmax, m_ctmin,
                            m_amct, 0, rans[1], rans[2]);
    wt *= CE.TChannelWeight(p[0] - p[5], p[1] - p[2], p[4], p[3], 0., m_alpha, 1., -1., m_amct, 0,
                            rans[3], rans[4]);
    if(wt != 0.) wt = 1.0 / wt / pow(2. * M_PI, 3 * 4. - 4.);
    // Mapper<Vec4D>::Print(__PRETTY_FUNCTION__, p, rans);
    return 1 / wt;
}

void C3_3::GeneratePoint(std::vector<Vec4D> &p, const std::vector<double> &ran) {
    GenerateNuclearPoint(p, ran);
    double s34_max = sqr((p[0] + p[1] - p[5]).Mass() - sqrt(s2));
    // TODO: set minimum from cuts
    // double s34_min = cuts->GetscutAmegic(std::string("34"));
    double s34_min = std::max(pow(sqrt(s3) + sqrt(s4), 2), 1E-8);
    Vec4D p34;
    double s34 = CE.MasslessPropMomenta(m_salpha, s34_min, s34_max, ran[0]);
    double tmass201 = Flavour((kf_code)(23)).Mass();
    CE.TChannelMomenta(p[0] - p[5], p[1], p[2], p34, s2, s34, tmass201, m_alpha, m_ctmax, m_ctmin,
                       m_amct, 0, ran[1], ran[2]);
    CE.TChannelMomenta(p[0] - p[5], p[1] - p[2], p[4], p[3], s4, s3, 0., m_alpha, 1., -1., m_amct,
                       0, ran[3], ran[4]);
    // Mapper<Vec4D>::Print(__PRETTY_FUNCTION__, p, ran);
}

double C3_3::GenerateWeight(const std::vector<Vec4D> &p, std::vector<double> &rans) {
    double wt = GenerateNuclearWeight(p, rans);
    double s34_max = sqr((p[0] + p[1] - p[5]).Mass() - sqrt(s2));
    // TODO: set minimum from cuts
    // double s34_min = cuts->GetscutAmegic(std::string("34"));
    double s34_min = std::max(pow(sqrt(s3) + sqrt(s4), 2), 1E-8);
    wt *= CE.MasslessPropWeight(m_salpha, s34_min, s34_max, dabs((p[3] + p[4]).Abs2()), rans[0]);
    double tmass201 = Flavour((kf_code)(23)).Mass();
    wt *= CE.TChannelWeight(p[0] - p[5], p[1], p[2], p[3] + p[4], tmass201, m_alpha, m_ctmax,
                            m_ctmin, m_amct, 0, rans[1], rans[2]);
    wt *= CE.TChannelWeight(p[0] - p[5], p[1] - p[2], p[4], p[3], 0., m_alpha, 1., -1., m_amct, 0,
                            rans[3], rans[4]);
    if(wt != 0.) wt = 1. / wt / pow(2. * M_PI, 3 * 4. - 4.);
    // Mapper<Vec4D>::Print(__PRETTY_FUNCTION__, p, rans);
    return 1 / wt;
}

void C3_4::GeneratePoint(std::vector<Vec4D> &p, const std::vector<double> &ran) {
    GenerateNuclearPoint(p, ran);
    double s34_max = sqr((p[0] + p[1] - p[5]).Mass() - sqrt(s2));
    // TODO: set minimum from cuts
    // double s34_min = cuts->GetscutAmegic(std::string("34"));
    double s34_min = std::max(pow(sqrt(s3) + sqrt(s4), 2), 1E-8);
    Vec4D p34;
    double s34 = CE.MasslessPropMomenta(m_salpha, s34_min, s34_max, ran[0]);
    CE.TChannelMomenta(p[0] - p[5], p[1], p[2], p34, s2, s34, 0., m_alpha, m_ctmax, m_ctmin, m_amct,
                       0, ran[1], ran[2]);
    CE.TChannelMomenta(p[0] - p[5], p[1] - p[2], p[3], p[4], s3, s4, 0., m_alpha, 1., -1., m_amct,
                       0, ran[3], ran[4]);
    // Mapper<Vec4D>::Print(__PRETTY_FUNCTION__, p, ran);
}

double C3_4::GenerateWeight(const std::vector<Vec4D> &p, std::vector<double> &rans) {
    double wt = GenerateNuclearWeight(p, rans);
    double s34_max = sqr((p[0] + p[1] - p[5]).Mass() - sqrt(s2));
    // TODO: set minimum from cuts
    // double s34_min = cuts->GetscutAmegic(std::string("34"));
    double s34_min = std::max(pow(sqrt(s3) + sqrt(s4), 2), 1E-8);
    wt *= CE.MasslessPropWeight(m_salpha, s34_min, s34_max, dabs((p[3] + p[4]).Abs2()), rans[0]);
    wt *= CE.TChannelWeight(p[0] - p[5], p[1], p[2], p[3] + p[4], 0., m_alpha, m_ctmax, m_ctmin,
                            m_amct, 0, rans[1], rans[2]);
    wt *= CE.TChannelWeight(p[0] - p[5], p[1] - p[2], p[3], p[4], 0., m_alpha, 1., -1., m_amct, 0,
                            rans[3], rans[4]);
    if(wt != 0.) wt = 1.0 / wt / pow(2. * M_PI, 3 * 4. - 4.);
    // Mapper<Vec4D>::Print(__PRETTY_FUNCTION__, p, rans);
    return 1 / wt;
}

void C3_5::GeneratePoint(std::vector<Vec4D> &p, const std::vector<double> &ran) {
    GenerateNuclearPoint(p, ran);
    double s34_max = sqr((p[0] + p[1] - p[5]).Mass() - sqrt(s2));
    // TODO: set minimum from cuts
    // double s34_min = cuts->GetscutAmegic(std::string("34"));
    double s34_min = std::max(pow(sqrt(s3) + sqrt(s4), 2), 1E-8);
    Vec4D p34;
    double s34 = CE.MasslessPropMomenta(m_salpha, s34_min, s34_max, ran[0]);
    double tmass201 = Flavour((kf_code)(23)).Mass();
    CE.TChannelMomenta(p[0] - p[5], p[1], p[2], p34, s2, s34, tmass201, m_alpha, m_ctmax, m_ctmin,
                       m_amct, 0, ran[1], ran[2]);
    CE.TChannelMomenta(p[0] - p[5], p[1] - p[2], p[3], p[4], s3, s4, 0., m_alpha, 1., -1., m_amct,
                       0, ran[3], ran[4]);
    // Mapper<Vec4D>::Print(__PRETTY_FUNCTION__, p, ran);
}

double C3_5::GenerateWeight(const std::vector<Vec4D> &p, std::vector<double> &rans) {
    double wt = GenerateNuclearWeight(p, rans);
    double s34_max = sqr((p[0] + p[1] - p[5]).Mass() - sqrt(s2));
    // TODO: set minimum from cuts
    // double s34_min = cuts->GetscutAmegic(std::string("34"));
    double s34_min = std::max(pow(sqrt(s3) + sqrt(s4), 2), 1E-8);
    wt *= CE.MasslessPropWeight(m_salpha, s34_min, s34_max, dabs((p[3] + p[4]).Abs2()), rans[0]);
    double tmass201 = Flavour((kf_code)(23)).Mass();
    wt *= CE.TChannelWeight(p[0] - p[5], p[1], p[2], p[3] + p[4], tmass201, m_alpha, m_ctmax,
                            m_ctmin, m_amct, 0, rans[1], rans[2]);
    wt *= CE.TChannelWeight(p[0] - p[5], p[1] - p[2], p[3], p[4], 0., m_alpha, 1., -1., m_amct, 0,
                            rans[3], rans[4]);
    if(wt != 0.) wt = 1.0 / wt / pow(2. * M_PI, 3 * 4. - 4.);
    // Mapper<Vec4D>::Print(__PRETTY_FUNCTION__, p, rans);
    return 1 / wt;
}

void C3_6::GeneratePoint(std::vector<Vec4D> &p, const std::vector<double> &ran) {
    GenerateNuclearPoint(p, ran);
    double s34_max = sqr((p[0] + p[1] - p[5]).Mass() - sqrt(s2));
    // TODO: set minimum from cuts
    // double s34_min = cuts->GetscutAmegic(std::string("34"));
    double s34_min = std::max(pow(sqrt(s3) + sqrt(s4), 2), 1E-8);
    Vec4D p34;
    double s34 = CE.MasslessPropMomenta(m_salpha, s34_min, s34_max, ran[0]);
    CE.TChannelMomenta(p[0] - p[5], p[1], p34, p[2], s34, s2, 0., m_alpha, m_ctmax, m_ctmin, m_amct,
                       0, ran[1], ran[2]);
    CE.Isotropic2Momenta(p34, s3, s4, p[3], p[4], ran[3], ran[4]);
    // Mapper<Vec4D>::Print(__PRETTY_FUNCTION__, p, ran);
}

double C3_6::GenerateWeight(const std::vector<Vec4D> &p, std::vector<double> &rans) {
    double wt = GenerateNuclearWeight(p, rans);
    double s34_max = sqr((p[0] + p[1] - p[5]).Mass() - sqrt(s2));
    // TODO: set minimum from cuts
    // double s34_min = cuts->GetscutAmegic(std::string("34"));
    double s34_min = std::max(pow(sqrt(s3) + sqrt(s4), 2), 1E-8);
    wt *= CE.MasslessPropWeight(m_salpha, s34_min, s34_max, dabs((p[3] + p[4]).Abs2()), rans[0]);
    wt *= CE.TChannelWeight(p[0] - p[5], p[1], p[3] + p[4], p[2], 0., m_alpha, m_ctmax, m_ctmin,
                            m_amct, 0, rans[1], rans[2]);
    wt *= CE.Isotropic2Weight(p[3], p[4], rans[3], rans[4], -1, 1);
    if(wt != 0.) wt = 1.0 / wt / pow(2. * M_PI, 3 * 4. - 4.);
    // Mapper<Vec4D>::Print(__PRETTY_FUNCTION__, p, rans);
    return 1 / wt;
}

void C3_7::GeneratePoint(std::vector<Vec4D> &p, const std::vector<double> &ran) {
    GenerateNuclearPoint(p, ran);
    double s34_max = sqr((p[0] + p[1] - p[5]).Mass() - sqrt(s2));
    // TODO: set minimum from cuts
    // double s34_min = cuts->GetscutAmegic(std::string("34"));
    double s34_min = std::max(pow(sqrt(s3) + sqrt(s4), 2), 1E-8);
    Flavour fl34 = Flavour((kf_code)(23));
    Vec4D p34;
    double s34 = CE.MassivePropMomenta(fl34.Mass(), fl34.Width(), 1, s34_min, s34_max, ran[0]);
    CE.TChannelMomenta(p[0] - p[5], p[1], p34, p[2], s34, s2, 0., m_alpha, m_ctmax, m_ctmin, m_amct,
                       0, ran[1], ran[2]);
    CE.Isotropic2Momenta(p34, s3, s4, p[3], p[4], ran[3], ran[4]);
    // Mapper<Vec4D>::Print(__PRETTY_FUNCTION__, p, ran);
}

double C3_7::GenerateWeight(const std::vector<Vec4D> &p, std::vector<double> &rans) {
    double wt = GenerateNuclearWeight(p, rans);
    double s34_max = sqr((p[0] + p[1] - p[5]).Mass() - sqrt(s2));
    // TODO: set minimum from cuts
    // double s34_min = cuts->GetscutAmegic(std::string("34"));
    double s34_min = std::max(pow(sqrt(s3) + sqrt(s4), 2), 1E-8);
    Flavour fl34 = Flavour((kf_code)(23));
    wt *= CE.MassivePropWeight(fl34.Mass(), fl34.Width(), 1, s34_min, s34_max,
                               dabs((p[3] + p[4]).Abs2()), rans[0]);
    wt *= CE.TChannelWeight(p[0] - p[5], p[1], p[3] + p[4], p[2], 0., m_alpha, m_ctmax, m_ctmin,
                            m_amct, 0, rans[1], rans[2]);
    wt *= CE.Isotropic2Weight(p[3], p[4], rans[3], rans[4], -1, 1);
    if(wt != 0.) wt = 1.0 / wt / pow(2. * M_PI, 3 * 4. - 4.);
    // Mapper<Vec4D>::Print(__PRETTY_FUNCTION__, p, rans);
    return 1 / wt;
}
