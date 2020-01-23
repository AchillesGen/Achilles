#ifndef UTILITIES_HH
#define UTILITIES_HH

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <vector>

const std::array<double, 3> ToCartesian(const std::array<double, 3>& vec);
bool sortPairSecond(const std::pair<std::size_t, double>& a,
                    const std::pair<std::size_t, double>& b);

class Brent {
    public:
        // Constructor
        Brent(std::function<double(double)> const& func, double tol=1e-6)
            : m_func(func), m_tol(tol) {}

        // Functions
        double CalcRoot(double, double) const;

    private:
        // Functions
        inline void swap(double &fa, double &fb, double &a, double &b) const {
            if(std::abs(fa) < std::abs(fb)) {
                double tmp = a;
                a = b;
                b = tmp;
                tmp = fa;
                fa = fb;
                fb = tmp;
            }
            return;
        }

        // Variables
        std::function<double(double)> m_func;
        double m_tol;
};

template<typename T>
class LogspaceGen {
private:
    T curValue, base;

public:
    LogspaceGen(T first, T base) : curValue(first), base(base) {}

    T operator()() {
        T retval = curValue;
        curValue *= base;
        return retval;
    }
};

inline std::vector<double> Logspace(double start, double stop, int num = 50, double base = 10) {
    double realStart = pow(base, start);
    double realBase = pow(base, (stop-start)/(num-1));

    std::vector<double> retval;
    retval.reserve(num);
    std::generate_n(std::back_inserter(retval), num, LogspaceGen<double>(realStart,realBase));
    return retval;
} 

template<typename T>
class LinspaceGen {
private:
    T curValue, base;

public:
    LinspaceGen(T first, T base) : curValue(first), base(base) {}

    T operator()() {
        T retval = curValue;
        curValue += base;
        return retval;
    }
};

inline std::vector<double> Linspace(double start, double stop, int num = 50) {
    double step = (stop-start)/(num-1);

    std::vector<double> retval;
    retval.reserve(num);
    std::generate_n(std::back_inserter(retval), num, LinspaceGen<double>(start,step));
    return retval;
} 

#endif // end of include guard: UTILITIES_HH
