#ifndef THREEVECTOR_HH
#define THREEVECTOR_HH

#include <array>
#include <cmath>
#include <iosfwd>

class ThreeVector {
    public:
        // Constructors and Destructors
        ThreeVector() noexcept : vec({0, 0, 0}) {}
        ThreeVector(std::array<double, 3> p) noexcept : vec(p) {}
        ThreeVector(double pX, double pY, double pZ) noexcept : vec({pX, pY, pZ}) {}
        ThreeVector(const ThreeVector& other) noexcept : vec(other.vec) {}
        ThreeVector(const ThreeVector&& other) noexcept : vec(std::move(other.vec)) {}
        ~ThreeVector() = default;

        // Setters
        void SetXYZ(const std::array<double, 3> p) noexcept {
            vec = p;
        }
        void SetXYZ(const double& x, const double& y, const double& z) noexcept {
            vec = std::array<double, 3>{x, y, z};
        }
        void SetPxPyPz(const std::array<double, 3> p) noexcept {
            vec = p;
        }
        void SetPxPyPz(const double& pX, const double& pY, const double& pZ) noexcept {
            vec = std::array<double, 3>{pX, pY, pZ};
        }
        void SetX(const double& x) noexcept {vec[0] = x;}
        void SetY(const double& y) noexcept {vec[1] = y;}
        void SetZ(const double& z) noexcept {vec[2] = z;}
        void SetPx(const double& pX) noexcept {vec[0] = pX;}
        void SetPy(const double& pY) noexcept {vec[1] = pY;}
        void SetPz(const double& pZ) noexcept {vec[2] = pZ;}

        // Getters
        const std::array<double, 3>& Position() const noexcept {return vec;}
        const double& X() const noexcept {return vec[0];}
        const double& Y() const noexcept {return vec[1];}
        const double& Z() const noexcept {return vec[2];}
        const double& Px() const noexcept {return vec[0];}
        const double& Py() const noexcept {return vec[1];}
        const double& Pz() const noexcept {return vec[2];}
        const double Pt2() const noexcept {return pow(vec[0], 2) + pow(vec[1], 2);}
        const double Pt() const noexcept {return sqrt(Pt2());}
        const double P2() const noexcept {return (*this)*(*this);} 
        const double P() const noexcept {return sqrt(P2());}
        const double Magnitude2() const noexcept {return P2();}
        const double Magnitude() const noexcept {return P();}
        const double Theta() const noexcept;
        const double Phi() const noexcept;

        // Functions
        const double Dot(const ThreeVector& other) const noexcept {
            return (*this) * other;
        }
        const ThreeVector Cross(const ThreeVector&) const noexcept;
        const ThreeVector Unit() const noexcept;
        const std::string ToString() const noexcept;

        // Operator Overloads
        ThreeVector& operator=(const ThreeVector& other) noexcept {
            vec = other.vec;
            return *this;
        }
        ThreeVector& operator+=(const ThreeVector&) noexcept;
        ThreeVector& operator-=(const ThreeVector&) noexcept;
        ThreeVector& operator*=(const double&) noexcept;
        ThreeVector& operator/=(const double&);

        double operator*(const ThreeVector&) const noexcept;
        ThreeVector operator-() const noexcept;
        ThreeVector operator+() const noexcept;
        ThreeVector operator*(const double&) const noexcept;
        ThreeVector operator/(const double&) const;
        ThreeVector operator-(const ThreeVector&) const noexcept;
        ThreeVector operator+(const ThreeVector&) const noexcept;

        // Comparison Operators
        bool operator==(const ThreeVector&) const noexcept;
        bool operator!=(const ThreeVector& other) const noexcept {return !(*this == other);}

        // Access Operators
        double& operator[](const std::size_t& idx) {
            if(idx > 2) throw std::range_error("Max value is 2.");
            return vec[idx];
        }
        const double& operator[](const std::size_t& idx) const {
            if(idx > 2) throw std::range_error("Max value is 2.");
            return vec[idx];
        }

        // Stream Operators
        friend std::ostream& operator<<(std::ostream&, const ThreeVector&);
        friend std::istream& operator>>(std::istream&, ThreeVector&);

    private:
        std::array<double, 3> vec;
};

ThreeVector operator*(const double&, const ThreeVector&) noexcept;

#endif // end of include guard: THREEVECTOR_HH
