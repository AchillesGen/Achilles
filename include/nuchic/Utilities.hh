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

namespace nuchic {

const std::array<double, 3> ToCartesian(const std::array<double, 3>& vec);
bool sortPairSecond(const std::pair<std::size_t, double>& a,
                    const std::pair<std::size_t, double>& b);

/// A class implementing the Brent method for root finding. This method of root finding
/// tries to take the optimal approach to finding the roots. Details on how the algorithm
/// works can be found at: https://en.wikipedia.org/wiki/Brent's_method
class Brent {
    public:
        /// Constructor
        ///@param func: Function to find the root of
        ///@param tol: The tolerance in the solution
        Brent(std::function<double(double)>  func, double tol=base_tol)
            : m_func(std::move(func)), m_tol(tol) {}

        // Functions
        /// Calculate the root in the given range. Throws a domain_error if no root exists
        ///@param a: The lower bound to search for the root
        ///@param b: The upper bound to search for the root
        ///@return double: The root that was found within the given tolerance
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
        static constexpr double base_tol = 1e-6;
};

/// Template class to act as a helper for generating ranges in log space. This is done in 
/// a way to mimic the logspace function found within the numpy python library.
template<typename T>
class LogspaceGen {
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
inline std::vector<double> Logspace(double start, double stop, std::size_t num = 50, double base = 10) {
    double realStart = pow(base, start);
    double realBase = pow(base, (stop-start)/(static_cast<double>(num)-1));

    std::vector<double> retval;
    retval.reserve(num);
    std::generate_n(std::back_inserter(retval), num, LogspaceGen<double>(realStart,realBase));
    return retval;
} 

/// Template class to act as a helper for generating ranges in linear space. This is done in 
/// a way to mimic the linspace function found within the numpy python library.
template<typename T>
class LinspaceGen {
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
    double step = (stop-start)/(static_cast<double>(num)-1);

    std::vector<double> retval;
    retval.reserve(num);
    std::generate_n(std::back_inserter(retval), num, LinspaceGen<double>(start,step));
    return retval;
} 

template<typename T>
constexpr T ipow(const T &x, const size_t &pow) {
    return pow == 0 ? 1 : x*ipow(x, pow-1);
}

template<typename T>
constexpr T factorial(const T &x) {
    return x == 0 ? 1 : x * factorial(x-1);
}

template<class ContainerT>
void tokenize(const std::string& str, ContainerT& tokens,
              const std::string& delimiters=" \n\t", bool trimEmpty=true) {
    std::string::size_type lastPos = 0, length = str.length();

    using value_type = typename ContainerT::value_type;
    using size_type = typename ContainerT::size_type;

    while(lastPos < length+1) {
        std::string::size_type pos = str.find_first_of(delimiters, lastPos);
        if(pos == std::string::npos) pos = length;

        if(pos != lastPos || !trimEmpty) {
            tokens.push_back(value_type(str.data()+lastPos,
                             static_cast<size_type>(pos)-lastPos));
        }
        lastPos = pos+1;
    }
}

namespace detail {

static const long Npts = 50;
const double bernoulli[Npts] = {
    1, -0.5               , 0.16666666666666666     , 0,
    -0.03333333333333333  , 0, 0.023809523809523808 , 0,
    -0.03333333333333333  , 0, 0.07575757575757576  , 0,
    -0.2531135531135531   , 0, 1.1666666666666667   , 0,
    -7.092156862745098    , 0, 54.971177944862156   , 0,
    -529.1242424242424    , 0, 6192.123188405797    , 0,
    -86580.25311355312    , 0, 1.4255171666666667e6 , 0,
    -2.7298231067816094e7 , 0, 6.015808739006424e8  , 0,
    -1.5116315767092157e10, 0, 4.296146430611667e11 , 0,
    -1.3711655205088334e13, 0, 4.883323189735932e14 , 0,
    -1.9296579341940068e16, 0, 8.416930475736827e17 , 0,
    -4.0338071854059454e19, 0, 2.1150748638081993e21, 0,
    -1.2086626522296526e23, 0  
};
const double fac_inv[Npts] = {
    1                     , 0.5                   , 0.16666666666666666   ,
    0.041666666666666664  , 0.008333333333333333  , 0.001388888888888889  ,
    0.0001984126984126984 , 0.0000248015873015873 , 2.7557319223985893e-6 ,
    2.755731922398589e-7  , 2.505210838544172e-8  , 2.08767569878681e-9   ,
    1.6059043836821613e-10, 1.1470745597729725e-11, 7.647163731819816e-13 ,
    4.779477332387385e-14 , 2.8114572543455206e-15, 1.5619206968586225e-16,
    8.22063524662433e-18  , 4.110317623312165e-19 , 1.9572941063391263e-20,
    8.896791392450574e-22 , 3.8681701706306835e-23, 1.6117375710961184e-24,
    6.446950284384474e-26 , 2.4795962632247972e-27, 9.183689863795546e-29 ,
    3.279889237069838e-30 , 1.1309962886447718e-31, 3.7699876288159054e-33,
    1.2161250415535181e-34, 3.800390754854744e-36 , 1.151633562077195e-37 ,
    3.387157535521162e-39 , 9.67759295863189e-41  , 2.688220266286636e-42 ,
    7.265460179153071e-44 , 1.911963205040282e-45 , 4.902469756513544e-47 ,
    1.225617439128386e-48 , 2.9893108271424046e-50, 7.117406731291439e-52 ,
    1.6552108677421951e-53, 3.7618428812322616e-55, 8.359650847182804e-57 ,
    1.817315401561479e-58 , 3.866628513960594e-60 , 8.055476070751238e-62 ,
    1.643974708316579e-63 , 3.2879494166331584e-65
};

}

double zeta(double s);
double eta(size_t n);
double PolyLog(long n, double x);

}

#endif // end of include guard: UTILITIES_HH
