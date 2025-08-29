#ifndef UTILITIES_HH
#define UTILITIES_HH

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace achilles {

// Bit operations taken from: https://graphics.stanford.edu/~seander/bithacks.html
inline unsigned int NextPermutation(unsigned int inp) {
    unsigned int t = inp | (inp - 1);
    return (t + 1) | (((~t & -~t) - 1) >> (__builtin_ctz(inp) + 1));
}

inline bool SetBit(unsigned int inp, unsigned int pos) {
    return inp & (1u << pos);
}

inline std::vector<unsigned int> SetBits(unsigned int inp, unsigned int size) {
    std::vector<unsigned int> set;
    for(unsigned int i = 0; i < size; ++i) {
        if(SetBit(inp, i)) set.push_back(inp & (1u << i));
    }
    return set;
}

inline bool IsPower2(unsigned int val) {
    return (val & (val - 1)) == 0;
}

inline unsigned int Log2(unsigned int val) {
    static const std::array<unsigned int, 5> b = {0xAAAAAAAA, 0xCCCCCCCC, 0xF0F0F0F0, 0xFF00FF00,
                                                  0xFFFF0000};
    unsigned int r = (val & b[0]) != 0;
    r |= static_cast<unsigned int>((val & b[4]) != 0) << 4;
    r |= static_cast<unsigned int>((val & b[3]) != 0) << 3;
    r |= static_cast<unsigned int>((val & b[2]) != 0) << 2;
    r |= static_cast<unsigned int>((val & b[1]) != 0) << 1;
    return r;
}

template <class ContainerT>
void tokenize(const std::string_view &str, ContainerT &tokens,
              const std::string_view &delimiters = " \n\t", bool trimEmpty = true) {
    std::string::size_type lastPos = 0, length = str.length();

    using value_type = typename ContainerT::value_type;
    using size_type = typename ContainerT::size_type;

    while(lastPos < length + 1) {
        std::string::size_type pos = str.find_first_of(delimiters, lastPos);
        if(pos == std::string::npos) pos = length;

        if(pos != lastPos || !trimEmpty) {
            tokens.push_back(
                value_type(str.data() + lastPos, static_cast<size_type>(pos) - lastPos));
        }
        lastPos = pos + 1;
    }
}

double ParseFraction(const std::string &);

constexpr int LeviCivita(const int i, const int j, const int k, const int l) {
    return (i == j || i == k || i == l || j == k || j == l || k == l)
               ? 0
               : (i - j) * (i - k) * (i - l) * (j - k) * (j - l) * (k - l) / 12;
}

class FourVector;
bool CheckMasses(const std::vector<achilles::FourVector> &, const std::vector<double> &,
                 double = 1e-8);

const std::array<double, 3> ToCartesian(const std::array<double, 3> &vec);
bool sortPairSecond(const std::pair<std::size_t, double> &a,
                    const std::pair<std::size_t, double> &b);

/// A class implementing the Brent method for root finding. This method of root finding
/// tries to take the optimal approach to finding the roots. Details on how the algorithm
/// works can be found at: https://en.wikipedia.org/wiki/Brent's_method
class Brent {
  public:
    /// Constructor
    ///@param func: Function to find the root of
    ///@param tol: The tolerance in the solution
    Brent(std::function<double(double)> func, double tol = base_tol)
        : m_func(std::move(func)), m_tol(tol) {}

    // Functions
    /// Calculate the root in the given range. Throws a domain_error if no root exists
    ///@param a: The lower bound to search for the root
    ///@param b: The upper bound to search for the root
    ///@return double: The root that was found within the given tolerance
    double CalcRoot(double, double) const;

    /// Calculate the minimum in the given range.
    ///@param ax: The lower bound to search for the minimum
    ///@param bx: The upper bound to search for the minimum
    ///@param cx: The initial guess for the minimum, (default: will find a guess)
    ///@return double: The minimum that was found within the given tolerance
    double Minimize(const double &ax, const double &bx, const double &cx = 1E99) const;

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
    inline void shift(double &a, double &b, double &c, const double &d) const {
        a = b;
        b = c;
        c = d;
    }
    void GetBracket(double &a, double &b, double &c, double &fa, double &fb, double &fc) const;

    // Variables
    std::function<double(double)> m_func;
    double m_tol;
    static constexpr double base_tol = 1e-6;
    static const size_t itmax = 100;
    static constexpr double tiny = 1.0E-20;
    static constexpr double glimit = 100.0;
    static constexpr double gold = 1.618034;
    static constexpr double cgold = 0.3819660;
    static constexpr double eps = std::numeric_limits<double>::epsilon() * 1.0E-3;
};

/// Template class to act as a helper for generating ranges in log space. This is done in
/// a way to mimic the logspace function found within the numpy python library.
template <typename T> class LogspaceGen {
  private:
    T curValue, base;

  public:
    /// Construct a helper object for generating equally spaced points in log space
    ///@param first: The value to start the range at
    ///@param base: The base of the exponent
    LogspaceGen(T first, T base_) : curValue(first), base(base_) {}

    /// Get the next value in the log space chain
    ///@return T: The next value in the chain
    T operator()() {
        T retval = curValue;
        curValue *= base;
        return retval;
    }
};

/// Function to generate a range of numbers equally spaced in log space from base^start to
/// base^stop. This function mimics the logspace function found within the numpy python library.
///@param start: The starting exponent to generate from
///@param stop: The ending exponent to generate to
///@param num: The number of points to generate within the range
///@param base: The base value to be used for the range
///@return std::vector<double>: A vector containing equally spaced points in log space
inline std::vector<double> Logspace(double start, double stop, std::size_t num = 50,
                                    double base = 10) {
    double realStart = pow(base, start);
    double realBase = pow(base, (stop - start) / (static_cast<double>(num) - 1));

    std::vector<double> retval;
    retval.reserve(num);
    std::generate_n(std::back_inserter(retval), num, LogspaceGen<double>(realStart, realBase));
    return retval;
}

/// Template class to act as a helper for generating ranges in linear space. This is done in
/// a way to mimic the linspace function found within the numpy python library.
template <typename T> class LinspaceGen {
  private:
    T curValue, step;

  public:
    /// Construct a helper object for generating equally spaced points in linear space
    ///@param first: The value to start the range at
    ///@param step: The value to add at each step
    LinspaceGen(T first, T step_) : curValue(first), step(step_) {}

    /// Get the next value in the linear space chain
    ///@return T: The next value in the chain
    T operator()() {
        T retval = curValue;
        curValue += step;
        return retval;
    }
};

/// Function to generate a range of numbers equally spaced in linear space from start to
/// stop. This function mimics the linspace function found within the numpy python library.
///@param start: The starting value to generate from
///@param stop: The ending value to generate to
///@param num: The number of points to generate within the range
///@return std::vector<double>: A vector containing equally spaced points in linear space
inline std::vector<double> Linspace(double start, double stop, std::size_t num = 50) {
    double step = (stop - start) / (static_cast<double>(num) - 1);

    std::vector<double> retval;
    retval.reserve(num);
    std::generate_n(std::back_inserter(retval), num, LinspaceGen<double>(start, step));
    return retval;
}

template <typename T> constexpr T ipow(const T &x, const size_t &pow) {
    return pow == 0 ? 1 : x * ipow(x, pow - 1);
}

template <typename T> constexpr T factorial(const T &x) {
    return x == 0 ? 1 : x * factorial(x - 1);
}

template <typename T> constexpr bool IsZero(const T &x, const T &tol) {
    return std::abs(x) < tol;
}

template <typename T> T Sign(const T &a) {
    return a > 0 ? 1 : (a < 0 ? -1 : 0);
}

// From arxiv:2201.01678
inline double Li2(double x) {
    static constexpr std::array<double, 6> p{0.9999999999999999502e+0, -2.6883926818565423430e+0,
                                             2.6477222699473109692e+0, -1.1538559607887416355e+0,
                                             2.0886077795020607837e-1, -1.0859777134152463084e-2};
    static constexpr std::array<double, 7> q{1.0000000000000000000e+0, -2.9383926818565635485e+0,
                                             3.2712093293018635389e+0, -1.7076702173954289421e+0,
                                             4.1596017228400603836e-1, -3.9801343754084482956e-2,
                                             8.2743668974466659035e-4};
    static constexpr double ZETA2 = M_PI * M_PI / 6;
    double result{};
    double sign = 1;

    if(x < -1) {
        result += log(1 - x) * (0.5 * log(1 - x) - log(-x)) - ZETA2;
        x = 1 / (1 - x);
    } else if(x == -1) {
        return -ZETA2;
    } else if(x < 0) {
        sign = -1;
        result += -1. / 2 * ipow(log(1 - x), 2);
        x = x / (x - 1);
    } else if(x == 0) {
        return 0;
    } else if(x < 0.5) {
        sign = 1;
    } else if(x == 0.5) {
        return ZETA2 / 2 - ipow(log(2), 2) / 2;
    } else if(x < 1) {
        result += ZETA2 - log(x) * log(1 - x);
        x = 1 - x;
        sign = -1;
    } else if(x == 1) {
        return ZETA2;
    } else if(x < 2) {
        result += ZETA2 - log(x) * (log(1 - 1 / x) + 1. / 2 * log(x));
        x = 1 - 1 / x;
    } else if(x == 2) {
        return M_PI * M_PI / 4;
    } else {
        sign = -1;
        result += ZETA2 * 2 - 1. / 2 * ipow(log(x), 2);
    }

    double x2 = x * x;
    double x4 = x2 * x2;
    double num = p[0] + p[1] * x + (p[2] + p[3] * x) * x2 + (p[4] + p[5] * x) * x4;
    double den = q[0] + q[1] * x + (q[2] + q[2] * x) * x2 + (q[4] + q[5] * x + q[6] * x2) * x4;
    return result + sign * x * num / den;
}

// TODO: Figure out a better place for these functions (maybe with an improved Settings class?)
size_t LevenshteinDistance(const std::string_view &src, const std::string_view &tgt);
std::string_view GetSuggestion(const std::vector<std::string> &opts, const std::string &query);

} // namespace achilles

#endif // end of include guard: UTILITIES_HH
