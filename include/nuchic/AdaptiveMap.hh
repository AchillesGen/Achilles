#ifndef ADAPTIVE_MAP_HH
#define ADAPTIVE_MAP_HH

#include <array>
#include <limits>
#include <mutex>
#include <string>
#include <vector>
#include <iosfwd>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "yaml-cpp/yaml.h"
#pragma GCC diagnostic pop

namespace nuchic {

enum class AdaptiveMapSplit {
    half,
    third,
    quarter
};

class AdaptiveMap {
    public:
        AdaptiveMap() = default;
        AdaptiveMap(size_t dims, size_t bins) 
            : m_dims(std::move(dims)), m_bins(std::move(bins)) {
                m_hist.resize(dims * (bins+1));
                for(size_t i = 0; i < m_bins + 1; ++i) {
                    m_hist[i] = static_cast<double>(i)/static_cast<double>(bins);
                }

                for(size_t i = 1; i < m_dims; ++i) {
                    std::copy(m_hist.begin(), m_hist.begin() + static_cast<int>(bins + 1),
                              m_hist.begin() + static_cast<int>(i * (bins + 1)));
                }
            }

        // Serialization
        bool Deserialize(std::istream &in);
        bool Serialize(std::ostream &out) const;

        // Bin locations
        double lower_edge(size_t dim, size_t bin) const { return m_hist[dim*(m_bins+1) + bin]; }
        double upper_edge(size_t dim, size_t bin) const { return m_hist[dim*(m_bins+1) + bin + 1]; }
        double width(size_t dim, size_t bin) const { return upper_edge(dim, bin) - lower_edge(dim, bin); }
        size_t FindBin(size_t, double) const;

        // Map information
        std::vector<double> Edges(size_t dim) const { 
            return std::vector<double>(m_hist.begin() + static_cast<int>(dim*(m_bins+1)),
                                       m_hist.begin() + static_cast<int>((dim+1)*(m_bins+1)));
        }
        size_t Bins() const { return m_bins; }
        size_t Dims() const { return m_dims; }
        // Used for testing purposes
        std::vector<double> Hist() const { return m_hist; }
        std::vector<double>& Hist() { return m_hist; }

        // Generate random numbers
        double operator()(std::vector<double>&);
        double GenerateWeight(const std::vector<double>&) const;

        // Update histograms
        void Adapt(const double&, const std::vector<double>&);
        void Split(AdaptiveMapSplit split = AdaptiveMapSplit::half);
            

    private:
        std::vector<double> m_hist;
        size_t m_dims{}, m_bins{};
};

}

namespace YAML {

template<>
struct convert<nuchic::AdaptiveMap> {
    static Node encode(const nuchic::AdaptiveMap &rhs) {
        Node node;
        node["ndims"] = rhs.Dims();
        node["nbins"] = rhs.Bins();
        for(size_t i = 0; i < rhs.Dims(); ++i) {
            node["hists"].push_back(rhs.Edges(i));
            node["hists"][i].SetStyle(YAML::EmitterStyle::Flow);
        }
        return node;
    }

    static bool decode(const Node &node, nuchic::AdaptiveMap &rhs) {
        // Ensure node has entries for ndims, nbins, and edges
        if(node.size() != 3) return false;

        // Load dimensions and bins and initialize map
        auto ndims = node["ndims"].as<size_t>();
        auto nbins = node["nbins"].as<size_t>();
        rhs = nuchic::AdaptiveMap(ndims, nbins);
        std::vector<double> hist(ndims*(nbins+1));

        // Ensure that the number of histograms matches the number of dims 
        if(node["hists"].size() != ndims) return false;
        for(size_t i = 0; i < ndims; ++i) {
            // Check that the histogram has the correct number of bins
            if(node["hists"][i].size() != nbins+1) return false; 
            for(size_t j = 0; j <= nbins; ++j) {
                hist[i*(nbins+1) + j] = node["hists"][i][j].as<double>();
            }
        }

        // Load histogram into map
        rhs.Hist() = hist;
        
        return true;
    }
};

}

#endif
