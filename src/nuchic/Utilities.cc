#include "nuchic/Utilities.hh"
#include "fmt/format.h"

#include <stdexcept>
#include <map>

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
    double fa = m_func(a), fb = m_func(b), fc{};
    if(std::abs(fa) < m_tol) return a;
    if(std::abs(fb) < m_tol) return b;
    if(fa*fb >= 0) throw std::domain_error("No root in given range");
    swap(fa, fb, a, b);
    double c = a;
    bool m_flag = true;
    double s{}, d = 0;
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

// The implementation for the polylogarithm is inspired from: 
// https://github.com/Expander/polylogarithm/blob/master/src/Li.cpp
namespace nuchic {
namespace detail {

bool is_close(double a, double b, double eps = 1e-16) {
    return std::abs(a - b) < eps;
}

bool is_even(size_t n) { return n % 2 == 0; }

double binomial(size_t n, size_t k) {
    if(k > n) throw std::runtime_error("PolyLog: Binomial coefficient requries n < k");

    double result = 1.0;
    if(k > n-k) k = n-k;
    for(size_t i = 0; i < k; ++i) {
        result *= static_cast<double>(n - i);
        result /= static_cast<double>(i + 1);
    }

    return result;
}

// TODO: Speed this up by only calculating if p_i > max(p_0, ..., p_(i-1))
using TArray = std::array<double, Npts>;
TArray Xn(size_t p) {
    static std::vector<TArray> xn;
    if(xn.size() < p+1) xn.resize(p+1);
    else return xn[p];

    {
        TArray ar;
        for(size_t ni = 0; ni < Npts; ++ni) {
            ar[ni] = bernoulli[ni];
        }
        xn[0] = ar;
    }

    for(size_t pi = 1; pi <= p; ++pi) {
        TArray ar;
        for(size_t ni = 0; ni < Npts; ++ni) {
            double sum = 0.0;
            for(size_t k = 0; k <= ni; ++k) {
                sum += binomial(ni, k)*bernoulli[ni-k]/static_cast<double>(k+1)
                     * xn[pi-1][k];
            }
            ar[ni] = sum;
        }
        xn[pi] = ar;
    }

    return xn[p];
}

double PolyLog_negative(long n, double x) {
    if(is_close(x, 1)) return std::numeric_limits<double>::infinity();

    size_t order = static_cast<size_t>(-n);
    double result = 0.0;
    for(size_t k = 0; k <= order; ++k) {
        double sum = 0.0;
        for(size_t j = 0; j <= k; ++j) {
            const double sign = is_even(j) ? -1.0 : 1.0;
            sum += sign*binomial(k, j)*ipow(static_cast<double>(j+1), order);
        }
        result += ipow(-x/(1-x), k+1)*sum;
    }
    return result;
}

double PolyLog_naive_sum(size_t n, double x) {
    double sum = 0, sum_old;
    size_t k = 0;

    do {
        k++;
        sum_old = sum;
        sum += ipow(x, k)/ipow(static_cast<double>(k), n);
    } while(!is_close(sum, sum_old) && k < std::numeric_limits<size_t>::max() - 2);

    return sum;
}

double Harmonic(long n) {
    double sum = 0.0;
    for(long h = 1; h <= n; ++h) sum += 1.0/static_cast<double>(h);
    return sum;
}

double PolyLog_around_one(long n, double x) {
    const double mu = log(x);
    double sum = 0.0, sum_old = 0.0;
    long k = 0;

    do {
        if(k == n-1) {
            k++;
            continue;
        }
        sum_old = sum;
        sum += zeta(static_cast<double>(n-k))/factorial(static_cast<double>(k))*std::pow(mu, k);
        k++;
    } while(!is_close(sum, sum_old) && k < std::numeric_limits<long>::max() - 2);

    return std::pow(mu, n-1)/factorial(static_cast<double>(n-1))*(Harmonic(n-1)-log(-mu))+sum;
}

}
}

double nuchic::zeta(double s) {
    static std::map<double, double> cache;
    if(cache.find(s) != cache.end()) return cache[s];
    else {
#if __cpp_lib_math_special_functions >= 201603
        const double result = std::riemann_zeta(s);
        cache[s] = result;
        return result;
#else
        if(s == 1) return std::numeric_limits<double>::infinity();

        double sum = 0.0, sum_old = 0.0;
        size_t n = 0;

        if(s >= 12) {
            n = 1;
            do {
                sum_old = sum;
                sum += pow(static_cast<double>(n), -s);
                n++;
            } while(!detail::is_close(sum_old, sum) && n < std::numeric_limits<size_t>::max() - 2);
            cache[s] = sum;
            return sum;
        }

        do { 
            sum_old = sum;
            double sub_sum = 0.0;

            for(size_t k = 0; k <= n; ++k) {
                const auto sign = detail::is_even(k) ? 1.0 : -1.0;
                sub_sum += detail::binomial(n, k)*sign*std::pow(k+1, -s);
            }
            
            sum += sub_sum*std::pow(2, -static_cast<long>(n+1));
            n++;
        } while(!detail::is_close(sum_old, sum) && n < std::numeric_limits<size_t>::max() - 2);
        sum /= (1.0 - pow(2.0, 1-s));
        cache[s] = sum;
        return sum;
#endif
    }
}

double nuchic::eta(size_t n) {
    return (1.0 - ipow(2, 1-n))*nuchic::zeta(static_cast<double>(n));
}

double nuchic::PolyLog(long n, double x) {
    static std::map<std::pair<size_t, double>, double> cache;
    std::pair<size_t, double> point{n ,x};
    if(cache.find(point) != cache.end()) {
        return cache[point];
    } else {
        if(std::abs(x) - 1 > 1E-16) {
            std::string errmsg = fmt::format("PolyLog: x = {} out of range (-1, 1)", x);
            throw std::runtime_error(errmsg);
        }
        if(n < 0) {
            const double result = detail::PolyLog_negative(n, x);
            cache[point] = result;
            return result;
        } else if(n == 0) {
            if(detail::is_close(x, 1)) return std::numeric_limits<double>::infinity();
            return x/(1.0-x);
        } else if(n == 1) return -log(1-x);
        else if(detail::is_close(x, 0)) return 0;
        else if(detail::is_close(x, 1)) return zeta(static_cast<double>(n));
        else if(detail::is_close(x, -1)) return -eta(static_cast<size_t>(n));
        else if(n >= 12) return detail::PolyLog_naive_sum(static_cast<size_t>(n), x);
        else if(detail::is_close(x, 1, 1e-2)) {
            const double result = detail::PolyLog_around_one(n, x);
            cache[point] = result;
            return result;
        } else {
            double u = -log(1-x);
            double p = 1., sum = 0;
            const auto xn = detail::Xn(static_cast<size_t>(n-2));
            for(size_t k = 0; k < detail::Npts; ++k) {
                p *= u;
                sum += xn[k]*p*detail::fac_inv[k];
            }
            cache[point] = sum;
            return sum;
        }
    }
}
