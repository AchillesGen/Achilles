#ifndef ACHILLES_XSECSPLINE
#define ACHILLES_XSECSPLINE

#include "Achilles/Histogram.hh"
#include "Achilles/Interpolation.hh"

namespace achilles {

class XSecSpline {
  public:
    XSecSpline() = default;
    XSecSpline(const std::string &, double emin, double emax, size_t nbins);
    void Fill(double x, double wgt);
    // sigma(E) is the MC integral estimator restricted to a bin: the weight sum
    // divided by the total number of samples thrown (self-tracked as the Fill
    // count) times the bin width, i.e. height / Entries(). This is independent of
    // the sampling density, so it stays unbiased as VEGAS reshapes the energy
    // draws -- provided each spline is filled exactly once per MC sample (zeros
    // included), so its Fill count equals the sample count N.
    void InitializeSpline();
    double Interpolate(double) const;

    // Cache Splines
    bool SaveState(std::ostream &) const;
    bool LoadState(std::istream &);

  private:
    bool m_initialized{false};
    double m_emin{}, m_emax{};
    Histogram m_xsec; // summed MC weights per energy bin
    Interp1D m_interp;
};

using NuElmPair = std::pair<int, int>;

class GeometryXSecSpline {
  public:
    GeometryXSecSpline() = default;
    void AddSpline(int nu_pid, int elm_pid, const XSecSpline &spline);
    XSecSpline &GetSpline(int nu_pid, int elm_pid);
    bool Empty() const { return m_xsecs.empty(); }
    // Fill every contained spline once for a single MC sample. nu_weights holds
    // the summed weight for species that contributed at this phase-space point;
    // any species absent from the map is filled with zero so all splines see the
    // same number of samples.
    void FillSample(double energy, const std::map<int, double> &nu_weights);
    // Finalize every contained spline (each self-normalizes by its Fill count).
    void InitializeSplines();
    double TotalXSec(double energy, int nu_pid, std::vector<int> mat_pid) const;

    // Cache Splines
    bool SaveState(std::ostream &) const;
    bool LoadState(std::istream &);

  private:
    std::map<NuElmPair, XSecSpline> m_xsecs;
};

} // end namespace achilles

#endif // ACHILLES_XSECSPLINE
