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

using NuElmPair = std::pair<int, int>;

class GeometryXSecSpline {
  public:
    GeometryXSecSpline() = default;
    void AddSpline(int nu_pid, int elm_pid, const XSecSpline &spline);
    XSecSpline &GetSpline(int nu_pid, int elm_pid);
    double TotalXSec(double energy, int nu_pid, std::vector<int> mat_pid) const;

    // Cache Splines
    bool SaveState(std::ostream &) const;
    bool LoadState(std::istream &);

  private:
    std::map<NuElmPair, XSecSpline> m_xsecs;
};

} // end namespace achilles

#endif // ACHILLES_XSECSPLINE
