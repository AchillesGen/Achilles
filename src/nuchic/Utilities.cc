#include "nuchic/Utilities.hh"

const std::array<double, 3> nuchic::ToCartesian(const std::array<double, 3>& vec) {
    const double x = vec[0] * std::sin(vec[1]) * std::cos(vec[2]);
    const double y = vec[0] * std::sin(vec[1]) * std::sin(vec[2]);
    const double z = vec[0] * std::cos(vec[1]);

    return {x, y, z};
}

bool nuchic::sortPairSecond(const std::pair<std::size_t, double>& a,
                    const std::pair<std::size_t, double>& b) {

    return a.second < b.second;
}

// Use the Brent algorithm to calculate the root of a given function
double nuchic::Brent::CalcRoot(double a, double b) const {
    double fa = m_func(a), fb = m_func(b), fc;
    if(fa*fb >= 0) throw std::domain_error("No root in given range");
    swap(fa, fb, a, b);
    double c = a;
    bool m_flag = true;
    double s, d = 0;
    while(fb != 0 && std::abs(b - a) > m_tol) {
        fc = m_func(c);
        if(fa != fc && fb != fc) {
            s = (a*fb*fc)/((fa-fb)*(fa-fc))+(b*fa*fc)/((fb-fa)*(fb-fc))+(c*fa*fb)/((fc-fa)*(fc-fb));
        } else {
            s = b - fb*(b-a)/(fb-fa);
        }
        if((s > std::max((3*a+b)/4, b) || s < std::min((3*a+b)/4, b)) ||
                (m_flag && std::abs(s-b) >= std::abs(b-c)/2) ||
                (!m_flag && std::abs(s-b) >= std::abs(c-d)/2) ||
                (m_flag && std::abs(b-c) < m_tol) ||
                (!m_flag && std::abs(c-d) < m_tol)) {
            s = (a+b)/2;
            m_flag = true;
        } else {
            m_flag = false;
        }

        double fs = m_func(s);
        d = c;
        c = b;
        if(fa*fs < 0) {
            b = s;
            fb = fs;
        } else {
            a = s;
            fa = fs;
        }
        swap(fa, fb, a, b);
    }
    return b;
}
