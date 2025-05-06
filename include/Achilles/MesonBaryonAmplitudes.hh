#pragma once

#include "Achilles/Interpolation.hh"
#include "f_polynomial.hh"
#include <array>
#include <complex>
#include <map>
#include <vector>

struct PartialWave {
    size_t L, J, I, S;

    bool operator<(const PartialWave &rhs) const {
        if(L != rhs.L) return L < rhs.L;
        if(J != rhs.J) return J < rhs.J;
        if(I != rhs.I) return I < rhs.I;
        return S < rhs.S;
    }
};

struct PartialWaveAmp {
    achilles::Interp1D real;
    achilles::Interp1D imag;
};

struct PartialWaveAmpW {
    std::vector<double> invariant_mass;
    std::vector<PartialWaveAmp> pwa;
    double mass_spect, mass_min_res;

    double PMax(double) const;
};

class MBAmplitudes {
  public:
    MBAmplitudes();
    MBAmplitudes(const MBAmplitudes &) = default;

    void readANL(int, int);
    void readANLDelta();
    bool readANLDeltaBlock(std::ifstream &file);
    void printPWA(size_t iMB, size_t fMB, int L, int iJ, int iI);

    double GetCG(size_t ichan, int twoJ, int twoIm, int twoIb) const {
        return CGcof[ichan][(twoJ - 1) / 2][twoIm + 2][twoIb + 2];
    }
    double CSTotal(double, size_t, int, int);

    int MesonPID_Cchan(size_t ichan) const { return Meson_PID_Cchan[ichan]; }
    int BaryonPID_Cchan(size_t ichan) const { return Baryon_PID_Cchan[ichan]; }
    size_t NChargeChannels() const { return nCchan; }

    int GetCchannel(int pidm, int pidb) const {
        for(size_t i = 0; i < nCchan; i++) {
            if((pidm == Meson_PID_Cchan[i]) && (pidb == Baryon_PID_Cchan[i])) {
                return static_cast<int>(i);
            }
        }
        return -1;
    }

    // TODO move some to private
    void initIso();

    void PIDlist();

    void initCChannels();

    void SetOpenCChannels();

    void Iso_Cchan(size_t iCChan, size_t &iMBchan, int &twoI3m, int &twoI3b);

    int iCChannel(size_t MBchan, int twoI3m, int twoI3b);

    size_t NMesonBaryonChannels() { return nMBchan; }

    size_t NFSChannels(size_t iCchan) { return NOpenCChannels[iCchan]; }

    std::vector<std::pair<double, size_t>> GetAllCSW(size_t ichan, double W) const;

    void SetCrossSectionsW();

    void CalcCrossSectionW_grid(size_t, size_t, double *);

    void SetLegendrecoeffs();
    void SetH_G_coeffs();

    double GetCSW(size_t i_i, size_t i_f, double W) const;
    void printCSW(size_t i_i, size_t i_f);

    f_Polynomial Get_CSpoly_W(double W, size_t i_i, size_t i_f) const;

  private:
    std::tuple<size_t, size_t> FromPartialWave(const std::string &) const;
    void SetThresholds();

    static constexpr int nPW = 20;
    std::string PWnames[nPW] = {"S11", "S31", "P11", "P13", "P31", "P33", "D13",
                                "D15", "D33", "D35", "F15", "F17", "F35", "F37",
                                "G17", "G19", "G37", "G39", "H19", "H39"};

    int twoJ_vec[nPW] = {1, 1, 1, 3, 1, 3, 3, 5, 3, 5, 5, 7, 5, 7, 7, 9, 7, 9, 9, 9};
    int L_vec[nPW] = {0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5};
    int twoI_vec[nPW] = {1, 3, 1, 1, 3, 3, 1, 1, 3, 3, 1, 1, 3, 3, 1, 1, 3, 3, 1, 3};

    static constexpr size_t maxMBchan = 4;
    size_t nMBchan = 0;
    std::array<std::string, maxMBchan> MBchannels = {"piN", "eta N", "K Lambda", "K Sigma"};
    std::array<double, maxMBchan> thresholds = {1080., 1490., 1615.,
                                                1690.};                // Minima in the datafiles
    std::array<double, maxMBchan> Mass_m = {138.5, 548.0, 495., 495.}; // Masses in ANL code
    std::array<double, maxMBchan> Mass_b = {938.5, 938.5, 1115.7, 1193.0};

    std::array<int, maxMBchan> twoIm_MBchannels = {2, 0, 1, 1};
    std::array<int, maxMBchan> twoIb_MBchannels = {1, 1, 0, 2};

    static constexpr int Wmax = 250;
    int num_W[maxMBchan][maxMBchan] = {{0}};

    double w_vec[maxMBchan][maxMBchan][Wmax];

    static constexpr int Lmax = 5;
    static constexpr int twoJmax = 9;
    std::complex<double> A_LJI[maxMBchan][maxMBchan][Wmax][Lmax + 1][2][2] = {{{{{{0.}}}}}};

    double CGcof[maxMBchan][2][5][5] = {{{{0.}}}}; // (I_3^{m}, I_3^{b}; Im, Ib | I , I_3) =
                                                   // CGcof[ichan][I-1/2][2I_3^m + 2][2I_3^b + 2]

    // Information for 'physical' (charge state) channels
    // These are numbered as in initchannels make
    static constexpr size_t maxchannels = 16; // piN(6) + eta N (2) + K Lamba (2) + K Sigma (6);

    size_t nCchan = 0;

    // int cC0_MBchannels[maxMBchan] = {0, 6, 8, 10};
    std::array<size_t, maxchannels> MBchan_Cchan;
    std::array<int, maxchannels> twoIm_Cchan;
    std::array<int, maxchannels> twoIb_Cchan;

    std::array<int, maxchannels> twoI3m_Cchan;
    std::array<int, maxchannels> twoI3b_Cchan;

    std::array<int, maxchannels> Baryon_PID_Cchan;
    std::array<int, maxchannels> Meson_PID_Cchan;

    std::array<size_t, maxchannels> NOpenCChannels = {0};
    size_t OpenCChannels[maxchannels][maxchannels] = {{0}};

    double CrossSectionsW[maxchannels][maxchannels][Wmax] = {{{0.}}};

    std::array<std::array<double, 12>, 12>
        A_Legendre; // Expansion coefficients of legendre polynomials

    // For angular distributions, expansion coefficients of H and G as polynomials
    double H_poly[Lmax + 1][Lmax + 1][2 * Lmax + 1] = {
        {{0.}}}; // H = P_n(x) * P_m(x)  = \sum a_k x^k, the coefficients a_k are H_poly[n][m] saved

    double G_poly[Lmax + 1][Lmax + 1][2 * Lmax + 1] = {{{0.}}}; // G = (1-x**2) P'_n(x) P'_m(x)

    std::map<PartialWave, PartialWaveAmpW> piDelta_pwa;
};
