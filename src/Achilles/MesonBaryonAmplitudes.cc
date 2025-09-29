#include <fstream>
#include <iostream>

#include "Achilles/Legendre.hh"
#include "Achilles/MesonBaryonAmplitudes.hh"
#include "Achilles/ParticleInfo.hh"
#include "Achilles/Utilities.hh"
#include "Achilles/f_polynomial.hh"

using achilles::ParticleInfo;
using achilles::PID;

double PartialWaveAmpW::PMax(double w) const {
    return 1 / (2 * w) *
           std::sqrt((pow(w * w - pow(mass_min_res, 2) - mass_spect * mass_spect, 2) -
                      4 * pow(mass_min_res, 2) * mass_spect * mass_spect));
}

MBAmplitudes::MBAmplitudes() {
    // Set the thresholds using Achilles ParticleInfo
    SetThresholds();

    for(int ii = 0; ii < 4; ii++) {
        for(int iff = 0; iff < 4; iff++) { readANL(ii, iff); }
    }
    readANLDelta();
    nMBchan = 4;
    initIso();
    initCChannels();
    PIDlist();
    SetOpenCChannels();
    SetCrossSectionsW();
    SetH_G_coeffs();
}

void MBAmplitudes::SetThresholds() {
    // Channel 1: piN
    thresholds[0] = std::max(
        ParticleInfo(PID::pionp()).Mass() + ParticleInfo(PID::neutron()).Mass(), thresholds[0]);

    // Channel 2: etaN
    thresholds[1] = std::max(ParticleInfo(PID::eta()).Mass() + ParticleInfo(PID::neutron()).Mass(),
                             thresholds[1]);

    // Channel 3: K Lambda
    thresholds[2] = std::max(
        ParticleInfo(PID::kaon0()).Mass() + ParticleInfo(PID::lambda0()).Mass(), thresholds[2]);

    // Channel 4: K Sigma
    thresholds[3] = std::max(ParticleInfo(PID::kaon0()).Mass() + ParticleInfo(PID::sigmam()).Mass(),
                             thresholds[3]);
}

void MBAmplitudes::readANL(int i_i, int i_f) {
    std::complex<double> C_i = -1.;
    C_i = std::sqrt(C_i);
    std::ifstream ing;
    std::string fnm = "ANL_" + std::to_string(i_i) + "-" + std::to_string(i_f) + ".dat";
    ing.open("./data/MesonBaryonAmplitudes/ANL/" + fnm);

    // skip the first two lines
    ing.ignore(1000, '\n');
    ing.ignore(1000, '\n');
    int iw = 0;
    while(!ing.eof()) {
        ing >> w_vec[i_i][i_f][iw];

        for(int iPW = 0; iPW < nPW; iPW++) {
            double Ar, Ai;
            ing >> Ar >> Ai;
            int L = L_vec[iPW];
            int twoJ = twoJ_vec[iPW];
            int iJ = 1;                      // J = L + 1/2
            if(twoJ - 2 * L < 0) { iJ = 0; } // J = L - 1/2

            int twoI = twoI_vec[iPW];
            int iI = (twoI - 1) / 2;

            A_LJI[i_i][i_f][iw][L][iJ][iI] = Ar + C_i * Ai;
        }

        iw++;
    }
    num_W[i_i][i_f] = iw - 1;

    ing.close();
}

void MBAmplitudes::readANLDelta() {
    std::ifstream input_file;
    input_file.open("./data/MesonBaryonAmplitudes/ANL/pwa-piDelta-pin.dat");
    std::string line;

    // Skip the first six lines
    for(size_t i = 0; i < 6; i++) { std::getline(input_file, line); }

    // Read all amplitude blocks
    while(readANLDeltaBlock(input_file)) {
        // Do nothing
    }
}

bool MBAmplitudes::readANLDeltaBlock(std::ifstream &input_file) {
    // Parse partial wave information
    std::string line;
    std::getline(input_file, line);
    std::vector<std::string> tokens;
    achilles::tokenize(line, tokens, " ");
    std::string pwa = tokens[2];
    size_t L = std::stoull(tokens[5]);
    double S = achilles::ParseFraction(tokens[7]);
    size_t spin = static_cast<size_t>(2 * S);
    auto [I, J] = FromPartialWave(pwa);

    if(piDelta_pwa.find(PartialWave{L, J, I, spin}) == piDelta_pwa.end()) {
        piDelta_pwa[PartialWave{L, J, I, spin}] = {};
        // std::cout << "Reading ANL partial wave " << pwa << " L = " << L << " S = " << S
        //           << std::endl;
    }

    // Read the W value
    std::getline(input_file, line);
    tokens.clear();
    achilles::tokenize(line, tokens, " ");
    // double W = std::stod(tokens[1]);

    // Read the amplitudes for each value of p
    std::vector<std::pair<double, std::complex<double>>> amplitudes;
    std::getline(input_file, line);
    while(std::getline(input_file, line)) {
        if(line.empty()) {
            // double pmax = 1.0 / (2 * W) *
            //               sqrt(pow(W * W - pow(938.5 + 138.5, 2) - 138.5 * 138.5, 2) -
            //                    4 * pow(938.5 + 138.5, 2) * 138.5 * 138.5);
            return true;
        }
        tokens.clear();
        achilles::tokenize(line, tokens, " ");
        double p = std::stod(tokens[0]);
        double Ar = std::stod(tokens[1]);
        double Ai = std::stod(tokens[2]);
        std::complex<double> A{Ar, Ai};
        amplitudes.push_back({p, A});
    }

    // double pmax = 1.0 / (2 * W) *
    //               sqrt(pow(W * W - pow(938.5 + 138.5, 2) - 138.5 * 138.5, 2) -
    //                    4 * pow(938.5 + 138.5, 2) * 138.5 * 138.5);
    // std::cout << "pmax = " << pmax << std::endl;

    return false;
}

void MBAmplitudes::printPWA(size_t iMB, size_t fMB, int L, int iJ, int iI) {
    std::cout << " PWA for " << MBchannels[iMB] << " - > " << MBchannels[fMB] << std::endl;
    int twoJ = 2 * L - 1 + iJ * 2;
    std::cout << " L , twoJ, twoI " << L << "  " << twoJ << "  " << 2. * iI + 1 << std::endl;

    for(int iW = 0; iW < num_W[iMB][fMB]; iW++) {
        std::cout << w_vec[iMB][fMB][iW] << "  ";

        std::cout << A_LJI[iMB][fMB][iW][L][iJ][iI].real() << "  ";
        std::cout << A_LJI[iMB][fMB][iW][L][iJ][iI].imag() << std::endl;
    }
}

void MBAmplitudes::initIso() {
    // Indices are 2*I_3 + 2, i.e. -1 - > 0, -1/2 -> 1, 1/2 -> 3, 1 -> 4
    // Indices for total J = 1/2 or 3/2 = 0 or 1 = (J -1/2)
    // chan = 0 =  piN , (1, 1/2)
    int chan = 0;
    CGcof[chan][1][4][3] = 1.; // 1, 1/2
    CGcof[chan][0][4][3] = 0.; // 1, 1/2

    double CG13 = sqrt(1. / 3.);
    double CG23 = sqrt(2. / 3.);

    CGcof[chan][0][4][1] = CG23, CGcof[chan][1][4][1] = CG13; // 1, -1/2

    CGcof[chan][0][2][3] = -CG13, CGcof[chan][1][2][3] = CG23; // 0 , 1/2

    CGcof[chan][0][2][1] = CG13, CGcof[chan][1][2][1] = CG23; // 0. -1/2

    CGcof[chan][0][0][3] = -CG23, CGcof[chan][1][0][3] = CG13; // -1, 1/2

    CGcof[chan][1][0][1] = 1.; //-1, -1/2

    // chan = 1 = \eta N , (0 , 1/2)
    chan = 1;
    int twoim = 0;
    int twoib = -1;
    CGcof[chan][0][twoim + 2][twoib + 2] = 1.;
    twoib = 1;
    CGcof[chan][0][twoim + 2][twoib + 2] = 1.;

    // chan = 2 = K \Lambda = (1/2, 0)
    chan = 2;
    twoib = 0;
    twoim = -1;
    CGcof[chan][0][twoim + 2][twoib + 2] = 1.;
    twoim = 1;
    CGcof[chan][0][twoim + 2][twoib + 2] = 1.;

    // chan = 3, K \Siga = (1/2, 1) -> (-1)^{J - 1/2 - 1} = -1 (J= 1/2), 1 (J = 3/2)
    chan = 3;
    for(int itauK = 1; itauK <= 3; itauK += 2) {
        for(int itauS = 0; itauS <= 4; itauS += 2) {
            CGcof[chan][0][itauK][itauS] = -1. * CGcof[0][0][itauS][itauK];
            CGcof[chan][1][itauK][itauS] = CGcof[0][1][itauS][itauK];
        }
    }
}

double MBAmplitudes::CSTotal(double W, size_t ichan, int twoIm, int twoIb) {
    double CGs[2] = {0., 0.};
    CGs[0] = GetCG(ichan, 1, twoIm, twoIb);
    CGs[1] = GetCG(ichan, 3, twoIm, twoIb);
    CGs[0] = CGs[0] * CGs[0];
    CGs[1] = CGs[1] * CGs[1];
    if(CGs[0] < 1e-3 && CGs[1] < 1e-3) { return 0.; }

    if(W < thresholds[ichan]) { return 0.; }

    int nW = num_W[ichan][ichan];
    if(W > w_vec[ichan][ichan][nW]) { return 0.; }

    // For lineair interpolation on equally spaced grid:
    double W_step = w_vec[ichan][ichan][3] - w_vec[ichan][ichan][2];
    int iW_m = static_cast<int>(std::floor(W / W_step));
    double W_m = w_vec[ichan][ichan][iW_m];
    double cofmax = (W - W_m) / W_step;
    double cofmin = 1. - cofmax;

    double CStot = 0.;

    for(int iI = 0; iI < 2; iI++) {
        if(CGs[iI] < 1e-3) { continue; }

        for(int iL = 0; iL < Lmax + 1; iL++) {
            for(int iJ = 0; iJ < 2; iJ++) {
                double Jfac = 2. * iL + iJ * 2.; // 2J + 1
                CStot += CGs[iI] * Jfac * (A_LJI[ichan][ichan][iW_m][iL][iJ][iI]).imag() * cofmin;
                CStot +=
                    CGs[iI] * Jfac * (A_LJI[ichan][ichan][iW_m + 1][iL][iJ][iI]).imag() * cofmax;
            }
        }
    }

    double PF = -4.; // times other stuff, fix !
    return CStot * PF;
}

void MBAmplitudes::initCChannels() {
    size_t iCChan = 0;

    for(size_t iMB = 0; iMB < nMBchan; iMB += 1) {
        int twoI_m = twoIm_MBchannels[iMB];
        int twoI_b = twoIb_MBchannels[iMB];
        for(int twoI3_m = -twoI_m; twoI3_m <= twoI_m; twoI3_m += 2) {
            for(int twoI3_b = -twoI_b; twoI3_b <= twoI_b; twoI3_b += 2) {
                twoI3m_Cchan[iCChan] = twoI3_m;
                twoI3b_Cchan[iCChan] = twoI3_b;

                twoIm_Cchan[iCChan] = twoI_m;
                twoIb_Cchan[iCChan] = twoI_b;

                MBchan_Cchan[iCChan] = iMB;

                iCChan++;
            }
        }
    }

    nCchan = iCChan;
}

void MBAmplitudes::PIDlist() {
    // Make a list of mesons and baryons corresponding to the channels
    // This is hardcoded, based on the PIDs (which will likely not change) but could be taken from
    // the parameterlist in principe
    Meson_PID_Cchan[0] = -211;
    Meson_PID_Cchan[1] = -211;
    Meson_PID_Cchan[2] = 111;
    Meson_PID_Cchan[3] = 111;
    Meson_PID_Cchan[4] = 211;
    Meson_PID_Cchan[5] = 211;

    Meson_PID_Cchan[6] = 221;
    Meson_PID_Cchan[7] = 221;

    Meson_PID_Cchan[8] = 311; // K0, will this complicate life ?
    Meson_PID_Cchan[9] = 321;

    Meson_PID_Cchan[10] = 311;
    Meson_PID_Cchan[11] = 311;
    Meson_PID_Cchan[12] = 311;
    Meson_PID_Cchan[13] = 321;
    Meson_PID_Cchan[14] = 321;
    Meson_PID_Cchan[15] = 321;

    Baryon_PID_Cchan[0] = 2112;
    Baryon_PID_Cchan[1] = 2212;
    Baryon_PID_Cchan[2] = 2112;
    Baryon_PID_Cchan[3] = 2212;
    Baryon_PID_Cchan[4] = 2112;
    Baryon_PID_Cchan[5] = 2212;

    Baryon_PID_Cchan[6] = 2112;
    Baryon_PID_Cchan[7] = 2212;

    Baryon_PID_Cchan[8] = 3122;
    Baryon_PID_Cchan[9] = 3122;

    Baryon_PID_Cchan[10] = 3112;
    Baryon_PID_Cchan[11] = 3212;
    Baryon_PID_Cchan[12] = 3222;
    Baryon_PID_Cchan[13] = 3112;
    Baryon_PID_Cchan[14] = 3212;
    Baryon_PID_Cchan[15] = 3222;
}

int MBAmplitudes::iCChannel(size_t MBchan, int twoI3m, int twoI3b) {
    // Return the index of the physical charge channel based on the charge states
    int ic0 = 0;
    for(size_t iMB = 0; iMB < MBchan; iMB++) {
        ic0 += (twoIm_MBchannels[iMB] + 1) * (twoIb_MBchannels[iMB] + 1);
    }

    int twoIm = twoIm_MBchannels[MBchan];
    int twoIb = twoIb_MBchannels[MBchan];
    int Nb = twoIb + 1;
    int im = (twoI3m + twoIm) / 2;
    int ib = (twoI3b + twoIb) / 2;

    return ic0 + im * Nb + ib; // index of the charge channel
}

void MBAmplitudes::Iso_Cchan(size_t iCChan, size_t &iMBchan, int &twoI3m, int &twoI3b) {
    iMBchan = MBchan_Cchan[iCChan];
    twoI3m = twoI3m_Cchan[iCChan];
    twoI3b = twoI3b_Cchan[iCChan];
}

// void MBAmplitudes::Iso_Cchan(int iCChan, int &iMBchan, int &twoI3m, int &twoI3b)
//{
//	int c0 = 0;
//	int iMB = 0;
//	for (iMB = 0 ; iMB < nMBchan ; iMB++)
//	{
//		int nstates = (twoIm_MBchannels[iMB] + 1)*(twoIb_MBchannels + 1);
//		if (iCChan < c0 + nstates){break;}
//		c0+=nstates;
//	}
//
//	iMBchan=iMB;
//
//	int ic = iCChan - c0;
//	int Nb = twoIb_MBchannels[iMB]+1;
//	int ib = ic%Nb
//	int im = (ic - ib)/Nb;
//
//	twoI3m = im*2 - twoIm_MBchannels[iMB];
//	twoI3b = ib*2 - twoIm_MBchannels[iMB];
//
// }

void MBAmplitudes::SetOpenCChannels() {
    for(size_t iCChan = 0; iCChan < nCchan; iCChan++) {
        size_t iMB_i = MBchan_Cchan[iCChan];
        int twoI3m_i = twoI3m_Cchan[iCChan];
        int twoI3b_i = twoI3b_Cchan[iCChan];
        int twoI3_i = twoI3m_i + twoI3b_i;
        double CG1_i = GetCG(iMB_i, 1, twoI3m_i, twoI3b_i);
        double CG3_i = GetCG(iMB_i, 3, twoI3m_i, twoI3b_i);

        size_t Nopen = 0;
        for(size_t i_f = 0; i_f < nCchan; i_f++) {
            size_t iMB_f = MBchan_Cchan[i_f];
            int twoI3m_f = twoI3m_Cchan[i_f];
            int twoI3b_f = twoI3b_Cchan[i_f];
            if(twoI3_i != twoI3m_f + twoI3b_f) { continue; }

            double CG1_f = GetCG(iMB_f, 1, twoI3m_f, twoI3b_f);
            double CG3_f = GetCG(iMB_f, 3, twoI3m_f, twoI3b_f);
            if((std::abs(CG1_f * CG1_i) < 1.0 * 1e-4) && (std::abs(CG3_f * CG3_i) < 1.0 * 1e-4)) {
                continue;
            } // should not happen in principle

            // Channel is open, set information
            OpenCChannels[iCChan][Nopen] = i_f;
            Nopen++;
        }

        NOpenCChannels[iCChan] = Nopen;
    }
}

void MBAmplitudes::SetCrossSectionsW() {
    for(size_t ichan_i = 0; ichan_i < nCchan; ichan_i++) {
        size_t nOpenFS = NOpenCChannels[ichan_i];
        for(size_t iFS = 0; iFS < nOpenFS; iFS++) {
            size_t ichan_f = OpenCChannels[ichan_i][iFS];
            CalcCrossSectionW_grid(ichan_i, ichan_f, CrossSectionsW[ichan_i][ichan_f]);
        }
    }
}

void MBAmplitudes::CalcCrossSectionW_grid(size_t iIS, size_t iFS, double *CS) {
    size_t iMB_i = MBchan_Cchan[iIS];
    int twoI3m_i = twoI3m_Cchan[iIS];
    int twoI3b_i = twoI3b_Cchan[iIS];
    double CG1_i = GetCG(iMB_i, 1, twoI3m_i, twoI3b_i);
    double CG3_i = GetCG(iMB_i, 3, twoI3m_i, twoI3b_i);

    size_t iMB_f = MBchan_Cchan[iFS];
    int twoI3m_f = twoI3m_Cchan[iFS];
    int twoI3b_f = twoI3b_Cchan[iFS];
    double CG1 = GetCG(iMB_f, 1, twoI3m_f, twoI3b_f) * CG1_i;
    double CG3 = GetCG(iMB_f, 3, twoI3m_f, twoI3b_f) * CG3_i;

    double CG_I[2] = {CG1, CG3};

    int nW = num_W[iMB_i][iMB_f];

    int iI_min = 0;
    int iI_max = 2;
    //		if (std::abs(CG1) < 1e-4){iI_min = 1;}
    //		if (std::abs(CG3) < 1e-4){iI_max = 1;}

    std::complex<double> Amp_iso;

    double mM = Mass_m[iMB_i];
    double mB = Mass_b[iMB_i];

    double mM2 = mM * mM;
    double mB2 = mB * mB;

    double Mfac = 4. * mM2 * mB2;

    for(int iW = 0; iW < nW; iW++) {
        // Calculate W-dependent factors here or later:
        double W = w_vec[iMB_i][iMB_f][iW];
        double W2 = W * W;
        CS[iW] = 0.;

        double PF = std::pow(W * W - mM2 - mB2, 2) - Mfac;
        if(PF < 0) { continue; }

        for(int L = 0; L < Lmax + 1; L++) {
            for(int iJ = 0; iJ < 2; iJ++) {
                double Jfac = 2. * L + iJ * 2.; // 2J + 1

                if(Jfac - 1 > twoJmax) { continue; }

                Amp_iso = 0.;
                for(int iI = iI_min; iI < iI_max; iI++) {
                    Amp_iso += CG_I[iI] * A_LJI[iMB_i][iMB_f][iW][L][iJ][iI];
                }
                CS[iW] += Jfac * std::pow(std::abs(Amp_iso), 2);
            }
        }

        //			CS[iW] = 2.*M_PI * 4.*W2/PF * CS[iW]; // unit 1/MeV^2
        CS[iW] = 197.32 * 197.32 * 10 * 2. * M_PI * 4. * W2 / PF * CS[iW]; // unit mB
    }
}

double MBAmplitudes::GetCSW(size_t i_i, size_t i_f, double W) const {
    // Simple lineair interpolation on a fixed grid
    int twoI3_i = twoI3m_Cchan[i_i];
    twoI3_i += twoI3b_Cchan[i_i];

    int twoI3_f = twoI3m_Cchan[i_f];
    twoI3_f += twoI3b_Cchan[i_f];

    if(twoI3_i != twoI3_f) { return 0.; }

    size_t iMB_i = MBchan_Cchan[i_i];
    size_t iMB_f = MBchan_Cchan[i_f];

    if(W < thresholds[iMB_i] || W < thresholds[iMB_f]) { return 0.; }

    int nW = num_W[iMB_i][iMB_f];
    if(W > w_vec[iMB_i][iMB_f][nW - 1]) { return 0.; }

    // For lineair interpolation on equally spaced grid:
    double W_step = w_vec[iMB_i][iMB_f][3] - w_vec[iMB_i][iMB_f][2];
    int iW_m = static_cast<int>(std::floor((W - w_vec[iMB_i][iMB_f][0]) / W_step));
    double W_m = w_vec[iMB_i][iMB_f][iW_m];
    double cofmax = (W - W_m) / W_step;
    double cofmin = 1. - cofmax;

    double CS_W = CrossSectionsW[i_i][i_f][iW_m] * cofmin;
    CS_W += CrossSectionsW[i_i][i_f][iW_m + 1] * cofmax;
    return CS_W;
}

std::vector<std::pair<double, size_t>> MBAmplitudes::GetAllCSW(size_t ichan, double W) const {
    std::vector<std::pair<double, size_t>> CS_fchan;
    if(ichan >= maxchannels) { return CS_fchan; } // return empty }
    size_t Nfchan = NOpenCChannels[ichan];

    for(size_t ifchan = 0; ifchan < Nfchan; ifchan++) {
        size_t fchan = OpenCChannels[ichan][ifchan];
        double CS = GetCSW(ichan, fchan, W);
        CS_fchan.push_back({CS, fchan});
    }

    return CS_fchan;
}

void MBAmplitudes::printCSW(size_t i_i, size_t i_f) {
    size_t iMB_i = MBchan_Cchan[i_i];
    size_t iMB_f = MBchan_Cchan[i_f];

    int nW = num_W[iMB_i][iMB_f];

    for(int iW = 0; iW < nW; iW++) {
        std::cout << w_vec[iMB_i][iMB_f][iW] << "  " << CrossSectionsW[i_i][i_f][iW] << std::endl;
    }
}

void MBAmplitudes::SetLegendrecoeffs() {
    LegendreCoefficients(Lmax, A_Legendre);
}

void MBAmplitudes::SetH_G_coeffs() {
    SetLegendrecoeffs();

    for(size_t Ln = 0; Ln < static_cast<size_t>(Lmax + 1); Ln++) {
        // for (int Lm = Ln ; Lm < Lmax+1; Lm++) // in principle symmetric, but can just do all of
        // them since this function is only called once if used correctly
        for(size_t Lm = 0; Lm < static_cast<size_t>(Lmax + 1); Lm++) {
            H_G_Coefficients(Ln, Lm, A_Legendre, H_poly[Ln][Lm], G_poly[Ln][Lm]);
        }
    }
}

f_Polynomial MBAmplitudes::Get_CSpoly_W(double W, size_t i_i, size_t i_f) const {
    size_t iMB_i = MBchan_Cchan[i_i];
    size_t iMB_f = MBchan_Cchan[i_f];

    //	if (W < thresholds[iMB_i] || W < thresholds[iMB_f]){ return 0.;} //error handle later

    //	int nW = num_W[iMB_i][iMB_f];
    //	if (W > w_vec[iMB_i][iMB_f][nW-1]){ return 0.; }

    // NN interpolation on equally spaced grid:
    double W_step = w_vec[iMB_i][iMB_f][3] - w_vec[iMB_i][iMB_f][2];
    int iW_m = static_cast<int>(std::floor((W - w_vec[iMB_i][iMB_f][0]) / W_step));
    double W_m = w_vec[iMB_i][iMB_f][iW_m];
    double cofmax = (W - W_m) / W_step;
    if(cofmax > 0.5) { iW_m = iW_m + 1; }

    double A_poly[2 * Lmax + 1] = {0.};

    // isospin sum
    int twoI3m_i = twoI3m_Cchan[i_i];
    int twoI3b_i = twoI3b_Cchan[i_i];
    double CG1_i = GetCG(iMB_i, 1, twoI3m_i, twoI3b_i);
    double CG3_i = GetCG(iMB_i, 3, twoI3m_i, twoI3b_i);

    int twoI3m_f = twoI3m_Cchan[i_f];
    int twoI3b_f = twoI3b_Cchan[i_f];
    double CG1 = GetCG(iMB_f, 1, twoI3m_f, twoI3b_f) * CG1_i;
    double CG3 = GetCG(iMB_f, 3, twoI3m_f, twoI3b_f) * CG3_i;

    std::complex<double> A_LJ[Lmax + 1][2];
    for(int L = 0; L < Lmax + 1; L++) {
        for(int iJ = 0; iJ < 2; iJ++) {
            A_LJ[L][iJ] = A_LJI[iMB_i][iMB_f][iW_m][L][iJ][0] * CG1;
            A_LJ[L][iJ] += A_LJI[iMB_i][iMB_f][iW_m][L][iJ][1] * CG3;
        }
    }

    for(int Ln = 0; Ln < Lmax + 1; Ln++) {
        double Lnd = 1. * Ln;
        // diagonal element first: even polynomial
        double ampH = pow(std::abs(Lnd * A_LJ[Ln][0] + (Lnd + 1.) * A_LJ[Ln][1]), 2);
        double ampG = pow(std::abs(A_LJ[Ln][0] - A_LJ[Ln][1]), 2);

        for(int ik = 0; ik < 2 * Lmax + 1; ik++) {
            A_poly[ik] += ampH * H_poly[Ln][Ln][ik] + ampG * G_poly[Ln][Ln][ik];
        }

        // Off diagonal elements, summed over Ln Lm + LmLn terms
        for(int Lm = Ln + 1; Lm < Lmax + 1; Lm++) {
            double Lmd = 1. * Lm;
            double ampG2 =
                2. * std::real(std::conj(A_LJ[Ln][0] - A_LJ[Ln][1]) * (A_LJ[Lm][0] - A_LJ[Lm][1]));
            std::complex<double> AH = std::conj(Lnd * A_LJ[Ln][0] + (Lnd + 1.) * A_LJ[Ln][1]) *
                                      (Lmd * A_LJ[Lm][0] + (Lmd + 1.) * A_LJ[Lm][1]);
            double ampH2 = 2. * std::real(AH);

            for(int ik = (Lm + Ln) % 2; ik < Lm + Ln + 1; ik += 2) {
                A_poly[ik] += ampG2 * G_poly[Ln][Lm][ik] + ampH2 * H_poly[Ln][Lm][ik];
            }
        }
    }

    //	return f_Polynomial(A_poly, 2*Lmax+1);
    return f_Polynomial(A_poly, 2 * Lmax);
}

std::tuple<size_t, size_t> MBAmplitudes::FromPartialWave(const std::string &pwa) const {
    for(size_t i = 0; i < nPW; i++) {
        if(pwa == PWnames[i]) { return {twoI_vec[i], twoJ_vec[i]}; }
    }
    throw std::runtime_error("Unknown partial wave " + pwa);
}

double f_Polynomial::eval_f(double x) {
    // Horner method
    double fx = Coeffs[order];
    for(int k = order - 1; k > -1; k--) { fx = Coeffs[k] + fx * x; }

    return fx;
}

double f_Polynomial::eval_df(double x) {
    double dfx = order * Coeffs[order];
    for(int k = order - 1; k > 0; k--) { dfx = Coeffs[k] * k + dfx * x; }

    return dfx;
}

double f_Polynomial::eval_If(double x) {
    double Ifx = Coeffs[order] / (1. * order + 1.);
    for(int k = order - 1; k > -1; k--) { Ifx = Coeffs[k] / (1. * k + 1.) + Ifx * x; }

    return Ifx * x;
}

void f_Polynomial::SetCoeffs(const double *cof, int _order) {
    order = _order;
    for(int k = 0; k <= order; k++) { Coeffs[k] = cof[k]; }
}
