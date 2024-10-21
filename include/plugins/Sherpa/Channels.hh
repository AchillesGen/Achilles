#ifndef PLUGINS_CHANNELS_HH
#define PLUGINS_CHANNELS_HH

#include <functional>
#include <utility>

// #include "ATOOLS/Math/Vector.H"
#include "Achilles/Mapper.hh"
#include "Achilles/PhaseSpaceFactory.hh"

namespace COMIX {
class Single_Process;
}

namespace ATOOLS {
class Flavour;
}

namespace PHASIC {

class Channels : public achilles::Mapper<ATOOLS::Vec4D> {
  public:
    double GenerateWeight(const std::vector<ATOOLS::Vec4D> &, std::vector<double> &) override = 0;
    void GeneratePoint(std::vector<ATOOLS::Vec4D> &, const std::vector<double> &) override = 0;
    size_t NDims() const override = 0;
    static std::string Name() { return "Sherpa Final State"; }
};

struct ChannelNode {
    std::shared_ptr<ChannelNode> m_left = nullptr, m_right = nullptr;
    long int m_pid;
    unsigned int m_idx;
};

class GenChannel : public Channels,
                   achilles::RegistrablePS<Channels, GenChannel, std::vector<double>> {
  public:
    GenChannel(size_t npart, std::vector<double> s)
        : m_n{npart}, m_nout{npart - 2}, m_s{std::move(s)} {
        m_p.resize(1ul << m_n);
    }
    bool InitializeChannel(std::shared_ptr<ChannelNode>);
    void GeneratePoint(std::vector<ATOOLS::Vec4D> &p, const std::vector<double> &rans) override;
    double GenerateWeight(const std::vector<ATOOLS::Vec4D> &p, std::vector<double> &rans) override {
        iran = 0;
        for(size_t i = 0; i < m_n; ++i) { m_p[1ul << i] = p[i]; }
        FillMomenta(m_nodes.get());
        auto lid = static_cast<size_t>((1 << m_n) - 2);
        m_p[lid] = p[0];
        double wgt = BuildWeight(m_nodes.get(), lid, rans, 0);
        if(wgt != 0)
            wgt = 1.0 / wgt / pow(2.0 * M_PI, 3.0 * static_cast<double>(m_nout) - 4) *
                  pow(1000, static_cast<double>(2 * m_nout - 4));
        return 1.0 / wgt;
    }
    size_t NDims() const override { return 3 * m_nout - 4; }
    void WriteChannel() const;
    YAML::Node ToYAML() const override {
        YAML::Node result;
        result["Name"] = m_name;
        result["Masses"] = m_s;
        return result;
    }
    static std::unique_ptr<Channels> Construct(const std::vector<double> &s) {
        return std::make_unique<GenChannel>(s.size(), s);
    }

    std::string ToString() {
        auto lid = static_cast<size_t>((1 << m_n) - 2);
        auto result = PrintPoint(m_nodes.get(), lid, 0);
        return result.substr(0, result.size() - 1);
    }
    std::string PrintPoint(ChannelNode *cur, size_t lid, size_t depth) const;
    std::string PrintSPoint(ChannelNode *node, size_t id) const;

  private:
    constexpr static size_t m_rid = 2, m_lid = 1;
    constexpr static double m_salpha = 0.9;
    constexpr static double m_amct = 1.0, m_alpha = 0.9, m_ctmax = 1.0, m_ctmin = -1.0;
    ChannelNode *LocateNode(ChannelNode *node, size_t id);
    void BuildPoint(ChannelNode *node, size_t lid, const std::vector<double> &, size_t depth = 0);
    void BuildSPoint(ChannelNode *node, size_t id, const std::vector<double> &);
    double BuildWeight(ChannelNode *node, size_t lid, std::vector<double> &, size_t depth = 0);
    double BuildSWeight(ChannelNode *node, size_t id, std::vector<double> &);
    double SCut(size_t);
    size_t SId(size_t);
    double PropMomenta(ChannelNode *, size_t, double, double, double);
    double PropWeight(ChannelNode *, size_t, double, double, double, double &);
    void TChannelMomenta(ChannelNode *, size_t, size_t, size_t, const std::vector<double> &);
    void SChannelMomenta(ChannelNode *, size_t, size_t, size_t, const std::vector<double> &);
    double TChannelWeight(ChannelNode *, size_t, size_t, size_t, std::vector<double> &);
    double SChannelWeight(ChannelNode *, size_t, size_t, size_t, std::vector<double> &);
    void FillMomenta(ChannelNode *);

    size_t m_n, m_nout;
    std::vector<double> m_s;
    std::string m_name;
    std::shared_ptr<ChannelNode> m_nodes{};
    std::vector<ATOOLS::Vec4D> m_p;
    size_t iran{};
};

} // namespace PHASIC

#endif
