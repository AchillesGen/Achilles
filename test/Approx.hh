#ifndef APPROX_HH
#define APPROX_HH

#include <sstream>

#include "catch2/catch.hpp"

// Inspiration taken from: https://github.com/catchorg/Catch2/issues/1467#issuecomment-473928075 

template<typename T>
class IsVectorApprox : public Catch::MatcherBase<T> {
    public:
        IsVectorApprox(T t) : _vec(t) {}

        bool match(const T& t) const override {
            bool result = true;
            for(size_t i = 0; i < t.Size(); ++i) {
                result &= (t[i] == Approx(_vec[i])
                    .epsilon(_epsilon)
                    .margin(_margin));
            }
            return result;
        }

        IsVectorApprox& epsilon(double eps) {
            _epsilon = eps;
            return *this;
        }

        IsVectorApprox& margin(double eps) {
            _margin = eps;
            return *this;
        }

        std::string describe() const override {
            std::ostringstream oss;
            oss << "is approximately " << _vec;
            return oss.str();
        }

    private:
        T _vec;
        double _epsilon{std::numeric_limits<double>::epsilon() * 100};
        double _margin{0.0};
};

#endif
