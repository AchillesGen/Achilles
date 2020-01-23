#ifndef FOURVECTOR_HH
#define FOURVECTOR_HH

#include <array>
#include <cmath>
#include <iosfwd>

class ThreeVector;

class FourVector {
    public:
        // Constructors and Destructors
        FourVector() noexcept : vec({0, 0, 0, 0}) {}
        FourVector(std::array<double, 4> p) noexcept : vec(p) {}
        FourVector(double pX, double pY, double pZ, double E) noexcept : vec({pX, pY, pZ, E}) {}
        FourVector(const FourVector& other) noexcept : vec(other.vec) {}
        FourVector(const FourVector&& other) noexcept : vec(std::move(other.vec)) {}
        ~FourVector() = default;

        // Setters
        void SetPxPyPzE(const std::array<double, 4> p) noexcept {
            vec = p;
        }
        void SetPxPyPzE(const double& pX, const double& pY, const double& pZ, const double& E) noexcept {
            vec = std::array<double, 4>{pX, pY, pZ, E};
        }
        void SetVectM(const ThreeVector& p, const double& mass) noexcept;
        void SetPx(const double& pX) noexcept {vec[0] = pX;}
        void SetPy(const double& pY) noexcept {vec[1] = pY;}
        void SetPz(const double& pZ) noexcept {vec[2] = pZ;}
        void SetE(const double& E) noexcept {vec[3] = E;}

        // Getters
        const std::array<double, 4>& Momentum() const noexcept {return vec;}
        const double& Px() const noexcept {return vec[0];}
        const double& Py() const noexcept {return vec[1];}
        const double& Pz() const noexcept {return vec[2];}
        const double& E() const noexcept {return vec[3];}
        const double Pt2() const noexcept {return pow(vec[0], 2) + pow(vec[1], 2);}
        const double Pt() const noexcept {return sqrt(Pt2());}
        const double P2() const noexcept {return Pt2() + pow(vec[2], 2);}
        const double P() const noexcept {return sqrt(P2());}
        const double M2() const noexcept {return (*this)*(*this);} 
        const double M() const noexcept;
        const double Magnitude2() const noexcept {return M2();}
        const double Magnitude() const noexcept {return M();}
        const double Theta() const noexcept;
        const double Phi() const noexcept;
        const double Rapidity() const noexcept;
        const double DeltaR(const FourVector&) const noexcept;
        const ThreeVector Vec3() const noexcept;

        // Functions
        FourVector Boost(const ThreeVector&) noexcept;
        FourVector Boost(const double&, const double&, const double&) noexcept;
        const ThreeVector BoostVector() const noexcept;
        const double Dot(const FourVector& other) const noexcept {
            return (*this) * other;
        }
        const std::string ToString() const noexcept;

        // Operator Overloads
        FourVector& operator=(const FourVector& other) noexcept {
            vec = other.vec;
            return *this;
        }
        FourVector& operator+=(const FourVector&) noexcept;
        FourVector& operator-=(const FourVector&) noexcept;
        FourVector& operator*=(const double&) noexcept;
        FourVector& operator/=(const double&);

        double operator*(const FourVector&) const noexcept;
        FourVector operator-() const noexcept;
        FourVector operator+() const noexcept;
        FourVector operator*(const double&) const noexcept;
        FourVector operator/(const double&) const;
        FourVector operator-(const FourVector&) const noexcept;
        FourVector operator+(const FourVector&) const noexcept;

        // Comparison Operators
        bool operator==(const FourVector&) const noexcept;
        bool operator!=(const FourVector& other) const noexcept {return !(*this == other);}

        // Access Operators
        double& operator[](const std::size_t& idx) {return vec[idx];}
        const double& operator[](const std::size_t& idx) const {return vec[idx];}

        // Stream Operators
        friend std::ostream& operator<<(std::ostream&, const FourVector&);
        friend std::istream& operator>>(std::istream&, FourVector&);

    private:
        std::array<double, 4> vec;
};

FourVector operator*(const double&, const FourVector&) noexcept;

#endif /* end of include guard: FOURVECTOR_HH */
