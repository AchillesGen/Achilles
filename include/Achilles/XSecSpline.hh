#ifndef ACHILLES_XSECSPLINE
#define ACHILLES_XSECSPLINE

#include "Achilles/Histogram.hh"
#include "Achilles/Interpolation.hh"

namespace achilles {

class XSecSpline {
  public:
    XSecSpline(const std::string &, double emin, double emax, size_t nbins);
    void Fill(double x, double wgt);
    void InitializeSpline();
    double Interpolate(double) const;

    // Cache Splines
    bool SaveState(std::ostream &) const;
    bool LoadState(std::istream &);

  private:
    bool m_initialized{false};
    double m_emin, m_emax;
    Histogram m_xsec;
    Interp1D m_interp;
};

using NuElmPair = std::pair<size_t, size_t>;

class GeometryXSecSpline {
  public:
    GeometryXSecSpline() = default;
    void AddSpline(size_t nu_pid, size_t elm_pid, const XSecSpline &spline);
    double TotalXSec(double energy, size_t nu_pid, std::vector<size_t> mat_pid) const;

  private:
    std::map<NuElmPair, XSecSpline> m_xsecs;
};

} // end namespace achilles

#endif // ACHILLES_XSECSPLINE
