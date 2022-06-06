#ifndef APPROX_HH
#define APPROX_HH

#include <complex>
#include <sstream>

#include "Achilles/FourVector.hh"

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

class VectorComplexApprox : public Catch::MatcherBase<std::vector<std::complex<double>>> {
    public:
        VectorComplexApprox(std::vector<std::complex<double>> t) : _vec(t) {}

        bool match(const std::vector<std::complex<double>> &t) const override {
            bool result = true;
            for(size_t i = 0; i < t.size(); ++i) {
                result &= (t[i].real() == Approx(_vec[i].real())
                                                .epsilon(_epsilon)
                                                .margin(_margin));
            }
            return result;
        }

        VectorComplexApprox& epsilon(double eps) {
            _epsilon = eps;
            return *this;
        }

        VectorComplexApprox& margin(double margin) {
            _margin = margin;
            return *this;
        }

        std::string describe() const override {
            std::ostringstream oss;
            oss << "is approximately {";
            for(const auto &elm : _vec) {
                oss << elm.real() << (elm.imag() < 0 ? "-" : "+") << elm.imag() << ", ";
            }
            oss << "\b\b}";
            return oss.str();
        }

    private:
        std::vector<std::complex<double>> _vec;
        double _epsilon{std::numeric_limits<double>::epsilon() * 100};
        double _margin{0.0};
};

class FourVectorApprox : public Catch::MatcherBase<achilles::FourVector> {
    public:
        FourVectorApprox(achilles::FourVector t) : _vec(t) {}

        bool match(const achilles::FourVector &t) const override {
            bool result = true;
            for(size_t i = 0; i < 4; ++i) {
                result &= (t[i] == Approx(_vec[i]).epsilon(_epsilon).margin(_margin));
            }
            return result;
        }

        FourVectorApprox& epsilon(double eps) {
            _epsilon = eps;
            return *this;
        }

        FourVectorApprox& margin(double margin) {
            _margin = margin;
            return *this;
        }

        std::string describe() const override {
            std::ostringstream oss;
            oss << "is approximately {" << _vec << "}";
            return oss.str();
        }

    private:
        achilles::FourVector _vec;
        double _epsilon{std::numeric_limits<double>::epsilon() * 100};
        double _margin{0.0};
};

class AllFourVectorApprox : public Catch::MatcherBase<std::vector<achilles::FourVector>> {
    public:
        AllFourVectorApprox(std::vector<achilles::FourVector> t) : _vec(t) {}

        bool match(const std::vector<achilles::FourVector> &t) const override {
            bool result = true;
            for(size_t i = 0; i < t.size(); ++i) {
                FourVectorApprox vector_matcher(_vec[i]);
                result &= vector_matcher.epsilon(_epsilon).margin(_margin).match(t[i]);
            }
            return result;
        }

        AllFourVectorApprox& epsilon(double eps) {
            _epsilon = eps;
            return *this;
        }

        AllFourVectorApprox& margin(double margin) {
            _margin = margin;
            return *this;
        }

        std::string describe() const override {
            std::ostringstream oss;
            oss << "is approximately {";
            for(const auto &elm : _vec) {
                oss << elm << ", ";
            }
            oss << "\b\b}";
            return oss.str();
        }

    private:
        std::vector<achilles::FourVector> _vec;
        double _epsilon{std::numeric_limits<double>::epsilon() * 100};
        double _margin{0.0};
};

#endif
