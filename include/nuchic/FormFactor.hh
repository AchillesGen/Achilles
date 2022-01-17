#ifndef FORMFACTORY_HH
#define FORMFACTORY_HH

#include <array>
#include <memory>

#include "yaml-cpp/node/node.h"

namespace nuchic {

class FormFactor {
    public:
        struct Values {
            double Gep, Gen, Gmp, Gmn; 
            double F1s, F2s, F1v, F2v;
            double Ges, Gms, Gev, Gmv;
        };

        FormFactor() = default;
        FormFactor(const FormFactor&) = default;
        FormFactor(FormFactor&&) = default;
        FormFactor& operator=(const FormFactor&) = default;
        FormFactor& operator=(FormFactor&&) = default;
        virtual ~FormFactor() = default;

        virtual Values operator()(double) const = 0;

        static std::unique_ptr<FormFactor> Build(const YAML::Node&);

    protected:
        void Fill(double, Values*) const;
};

class Dipole : public FormFactor {
    public:
        Dipole(const YAML::Node&);

        Values operator()(double) const override;

    private:
        double lambda, MA, muP, muN;
};

class Kelly : public FormFactor {
    public:
        Kelly(const YAML::Node&);

        Values operator()(double) const override;

    private:
        double lambda, MA, muP, muN;
        std::array<double, 4> termsEp{};
        std::array<double, 2> termsEn{};
        std::array<double, 4> termsMp{};
        std::array<double, 4> termsMn{};

        double Parameterization(std::array<double, 4> A, double x) const {
            return (1 + A[0]*x)/(1 + A[1]*x + A[2]*x*x + A[3]*x*x*x); 
        }
};

class BBBA : public FormFactor {
    public:
        BBBA(const YAML::Node&);

        Values operator()(double) const override;

    private:
        double muP, muN;
        std::array<double, 4> numEp{}, denEp{};
        std::array<double, 4> numEn{}, denEn{};
        std::array<double, 4> numMp{}, denMp{};
        std::array<double, 4> numMn{}, denMn{};

        double Numerator(std::array<double, 4> A, double x) const {
            double result = 0;
            for(size_t i = A.size()-1; i > 0; --i) {
                result += A[i];
                result *= x;
            }

            return result + A[0];
        }

        double Denominator(std::array<double, 4> B, double x) const {
            double result = 0;
            for(size_t i = B.size(); i > 0; --i) {
                result += B[i-1];
                result *= x;
            }

            return result + 1;
        }
};

class ArringtonHill : public FormFactor {
    public:
        ArringtonHill(const YAML::Node&);

        Values operator()(double) const override;

    private:
        double muP, muN, tcut, t0;
        std::array<double, 13> epParams{}, enParams{}, mpParams{}, mnParams{};

        double ZExpand(std::array<double, 13> A, double z) const {
            double result = 0;
            for(size_t i = A.size()-1; i > 0; --i) {
                result += A[i];
                result *= z;
            }
            return result + A[0];
        }
};

}

#endif
