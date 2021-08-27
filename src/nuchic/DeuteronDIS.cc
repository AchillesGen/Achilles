#include "nuchic/DeuteronDIS.hh"
#include "nuchic/ParticleInfo.hh"

using nuchic::DeuteronDIS;

DeuteronDIS::DeuteronDIS(const std::string &spectral, const std::string& pdf, int iset)
    : m_spectral{spectral}, m_pdf{pdf, iset} {}

double DeuteronDIS::qvec(double x, double q2) const {
    return sqrt(q2+pow(omega(x, q2), 2));
}

double DeuteronDIS::qvec2(double x, double q2) const {
    return q2+pow(omega(x, q2), 2);
}

double DeuteronDIS::omega(double x, double q2) const {
    return q2/(mass*x);
}

double DeuteronDIS::omega_t(double x, double q2, double p) const {
    return omega(x, q2)-2*sqrt(p*p+mn2)+2*mn-0.002;
}

double DeuteronDIS::q2t(double x, double q2, double p) const {
    return qvec2(x, q2) - pow(omega_t(x, q2, p), 2);
}

double DeuteronDIS::phase_space(double p, bool dummy) const {
    if(dummy) return 1;
    return 2*M_PI*p*p*m_spectral(p)/pow(2*M_PI, 3)*mn/sqrt(p*p + mn2);
}

double DeuteronDIS::xn(double x, double q2, double p, double costheta) const {
    if(p == 0) return x;
    double ep = sqrt(p*p+mn2);
    double pdotq = ep*omega(x, q2) - p*costheta*qvec(x, q2);
    return q2/2/pdotq;
}

double DeuteronDIS::F2(double x, double q2, double p, double costheta) const {
    double sintheta2 = 1 - costheta*costheta;
    double ep = sqrt(p*p+mn2);
    double pdotq = ep*omega_t(x, q2, p) - p*costheta*qvec(x, q2);
    double x_n = xn(x, q2, p, costheta);

    if(x_n < 0 || x_n > 1) return 0;
    double f1p = m_pdf.F1(PID::proton(), x_n, q2t(x, q2, p));
    double f2p = m_pdf.F2(PID::proton(), x_n, q2t(x, q2, p));
    double f1n = m_pdf.F1(PID::neutron(), x_n, q2t(x, q2, p));
    double f2n = m_pdf.F2(PID::neutron(), x_n, q2t(x, q2, p));

    double term1 = omega(x, q2)/mn * q2/qvec2(x, q2)*(1-q2/q2t(x,q2,p));
    double term2 = q2*q2/pow(qvec2(x, q2)*mn, 2)*pow(ep+pdotq/q2t(x, q2, p)*omega_t(x, q2, p), 2);
    double term3 = q2/qvec2(x, q2)/2*p*p*sintheta2/mn2;

    double ps = phase_space(p);

    return ps*(term1*(f1p+f1n) + (term2+term3)*(f2p+f2n));
}

double DeuteronDIS::F1(double x, double q2, double p, double costheta) const {
    double sintheta2 = 1 - costheta*costheta;
    double ep = sqrt(p*p+mn2);
    double pdotq = ep*omega(x, q2) - p*costheta*qvec(x, q2);
    double x_n = q2/2/pdotq;

    if(x_n < 0 || x_n > 1) return 0;
    if(q2t(x, q2, p) < 0) return 0;
    double f1p = m_pdf.F1(PID::proton(), x_n, q2t(x, q2, p));
    double f2p = m_pdf.F2(PID::proton(), x_n, q2t(x, q2, p));
    double f1n = m_pdf.F1(PID::neutron(), x_n, q2t(x, q2, p));
    double f2n = m_pdf.F2(PID::neutron(), x_n, q2t(x, q2, p));

    double f_p = mass*(f1p/mn + f2p/(2*omega(x, q2))*p*p*sintheta2/mn2);
    double f_n = mass*(f1n/mn + f2n/(2*omega(x, q2))*p*p*sintheta2/mn2);

    double ps = phase_space(p);

    return ps*(f_p + f_n);
}

std::vector<double> DeuteronDIS::F2All(double x, double q2, double p, double costheta) const {
    double sintheta2 = 1 - costheta*costheta;
    double ep = sqrt(p*p+mn2);
    double pdotq = ep*omega_t(x, q2, p) - p*costheta*qvec(x, q2);
    double x_n = xn(x, q2, p, costheta);

    std::vector<double> results(m_pdf.NSets());
    if(x_n < 0 || x_n > 1) return results;
    auto f1p = m_pdf.F1All(PID::proton(), x_n, q2t(x, q2, p));
    auto f2p = m_pdf.F2All(PID::proton(), x_n, q2t(x, q2, p));
    auto f1n = m_pdf.F1All(PID::neutron(), x_n, q2t(x, q2, p));
    auto f2n = m_pdf.F2All(PID::neutron(), x_n, q2t(x, q2, p));

    double term1 = omega(x, q2)/mn * q2/qvec2(x, q2)*(1-q2/q2t(x,q2,p));
    double term2 = q2*q2/pow(qvec2(x, q2)*mn, 2)*pow(ep+pdotq/q2t(x, q2, p)*omega_t(x, q2, p), 2);
    double term3 = q2/qvec2(x, q2)/2*p*p*sintheta2/mn2;

    double ps = phase_space(p);

    for(size_t iset = 0; iset < results.size(); ++iset)
        results[iset] = ps*(term1*(f1p[iset]+f1n[iset]) + (term2+term3)*(f2p[iset]+f2n[iset]));
    return results;
}

std::vector<double> DeuteronDIS::F1All(double x, double q2, double p, double costheta) const {
    double sintheta2 = 1 - costheta*costheta;
    double ep = sqrt(p*p+mn2);
    double pdotq = ep*omega(x, q2) - p*costheta*qvec(x, q2);
    double x_n = q2/2/pdotq;

    std::vector<double> results(m_pdf.NSets());
    if(x_n < 0 || x_n > 1) return results;
    auto f1p = m_pdf.F1All(PID::proton(), x_n, q2t(x, q2, p));
    auto f2p = m_pdf.F2All(PID::proton(), x_n, q2t(x, q2, p));
    auto f1n = m_pdf.F1All(PID::neutron(), x_n, q2t(x, q2, p));
    auto f2n = m_pdf.F2All(PID::neutron(), x_n, q2t(x, q2, p));
    double ps = phase_space(p);

    for(size_t iset = 0; iset < results.size(); ++iset) {
        double f_p = mass*(f1p[iset]/mn + f2p[iset]/(2*omega(x, q2))*p*p*sintheta2/mn2);
        double f_n = mass*(f1n[iset]/mn + f2n[iset]/(2*omega(x, q2))*p*p*sintheta2/mn2);

        results[iset] = ps*(f_p + f_n);
    }

    return results;
}
