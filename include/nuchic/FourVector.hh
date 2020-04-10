#ifndef FOURVECTOR_HH
#define FOURVECTOR_HH

#include <array>
#include <cmath>
#include <iosfwd>

#include "spdlog/fmt/ostr.h"

namespace nuchic {

class ThreeVector;

/// @brief FourVector is a container for dealing with four vectors
///
/// The FourVector class provides an easy to use container to handle four
/// component vectors, such as four-position and four-momentum
class FourVector {
    public:
        /// @name Constructors and Destructors
        ///@{

        /// Create an empty FourVector object
        FourVector() noexcept : vec({0, 0, 0, 0}) {}
        /// Create a ThreeVector object with values given by p
        ///@param p: A std::array<double, 4> containing the values for the vector
        FourVector(std::array<double, 4> p) noexcept : vec(p) {}
        /// Create a FourVector object with values pX, pY, pZ, and E
        ///@param pX: The px value of the vector
        ///@param pY: The py value of the vector
        ///@param pZ: The pz value of the vector
        ///@param E: The E value of the vector
        FourVector(double pX, double pY, double pZ, double E) noexcept 
            : vec({pX, pY, pZ, E}) {}
        /// Create a FourVector object from a ThreeVector and an energy
        ///@param other: ThreeVector object containing 3-momentum information
        ///@param E: Energy for the given vector
        FourVector(const ThreeVector& other, const double& E) noexcept;
        /// Create a copy of a FourVector object
        ///@param other: The vector to be copied
        FourVector(const FourVector& other) noexcept = default;
        /// Move a FourVector object to another
        ///@param other: The vector to be moved
        FourVector(FourVector&&) noexcept = default;

        /// Default destructor
        ~FourVector() = default;
        ///@}

        /// @name Setters
        /// @{
        /// These functions provide access to setting the parameters of the FourVector

        /// Set momentum variable based on pased in array
        ///@param p: momentum to be stored
        void SetPxPyPzE(const std::array<double, 4> p) noexcept {
            vec = p;
        }

        /// Set momentum variable based on pased in array
        ///@param pX: momentum in the x-direction to be stored
        ///@param pY: momentum in the y-direction to be stored
        ///@param pZ: momentum in the z-direction to be stored
        ///@param E: energy to be stored
        void SetPxPyPzE(const double& pX, const double& pY, const double& pZ, const double& E) noexcept {
            vec = std::array<double, 4>{pX, pY, pZ, E};
        }

        /// Set the four momentum variable given a three momentum and a mass
        ///@param p: The three momentum to use
        ///@param mass: The invariant mass of the object
        void SetVectM(const ThreeVector& p, const double& mass) noexcept;

        /// Set only the x momentum
        ///@param pX: momentum in the x-direction to be stored
        void SetPx(const double& pX) noexcept {vec[0] = pX;}

        /// Set only the y momentum
        ///@param pY: momentum in the y-direction to be stored
        void SetPy(const double& pY) noexcept {vec[1] = pY;}

        /// Set only the z momentum
        ///@param pZ: momentum in the z-direction to be stored
        void SetPz(const double& pZ) noexcept {vec[2] = pZ;}

        /// Set only the energy
        ///@param E: the energy
        void SetE(const double& E) noexcept {vec[3] = E;}
        ///@}

        /// @name Getters
        /// @{
        /// These functions provide get specific features from the FourVector object

        /// Return the momentum as an array
        ///@return std::array<double, 4>: An array containing the four momentum
        const std::array<double, 4>& Momentum() const noexcept {return vec;}

        /// Return the x-coordinate
        ///@return double: the x-coordinate
        const double& X() const noexcept {return vec[0];}

        /// Return the y-coordinate
        ///@return double: the y-coordinate
        const double& Y() const noexcept {return vec[1];}

        /// Return the z-coordinate
        ///@return double: the z-coordinate
        const double& Z() const noexcept {return vec[2];}

        /// Return the time
        ///@return double: the time
        const double& T() const noexcept {return vec[3];}

        /// Return the momentum in the x-direction
        ///@return double: the momentum in the x-direction
        const double& Px() const noexcept {return vec[0];}

        /// Return the momentum in the y-direction
        ///@return double: the momentum in the y-direction
        const double& Py() const noexcept {return vec[1];}

        /// Return the momentum in the z-direction
        ///@return double: the momentum in the z-direction
        const double& Pz() const noexcept {return vec[2];}

        /// Return the energy
        ///@return double: the energy
        const double& E() const noexcept {return vec[3];}

        /// Return the transverse momentum squared
        ///@return double: Transverse momentum squared
        double Pt2() const noexcept {return pow(vec[0], 2) + pow(vec[1], 2);}

        /// Return the transverse momentum
        ///@return double: Transverse momentum
        double Pt() const noexcept {return sqrt(Pt2());}

        /// Return the three momentum squared
        ///@return double: Three momentum squared
        double P2() const noexcept {return Pt2() + pow(vec[2], 2);}

        /// Return the three momentum
        ///@return double: Three momentum
        double P() const noexcept {return sqrt(P2());}

        /// Return the invariant mass squared
        ///@return double: The invariant mass squared
        double M2() const noexcept {return (*this)*(*this);} 

        /// Return the invariant mass
        ///@return double: The invariant mass
        double M() const noexcept;

        /// Return the Minkowski magnitude squared
        ///@return double: The magnitude squared
        double Magnitude2() const noexcept {return M2();}

        /// Return the Minkowski magnitude
        ///@return double: The magnitude
        double Magnitude() const noexcept {return M();}

        /// Return the angle between the z-axis and the transverse plane
        ///@return double: Angle between z-axis and transverse plane
        double Theta() const noexcept;

        /// Return the angle in the transverse plane
        ///@return double: The angle in the transverse plane
        double Phi() const noexcept;

        /// Return the rapidity of the momentum
        ///@return double: The rapidity associated with the momentum
        double Rapidity() const noexcept;

        /// Return the distance in the Eta-Phi plane between two four vectors
        ///@return double: Distance between two four vectors in Eta-Phi plane
        double DeltaR(const FourVector&) const noexcept;

        /// Return the three momentum component
        ///@return ThreeVector: The three momentum component
        ThreeVector Vec3() const noexcept;
        ///@}

        /// @name Functions
        /// @{

        /// Boost the four vector to the frame given by the three vector velocities
        ///@param beta: The boost vector to determine the frame
        ///@return FourVector: The vector in the corresponding frame
        FourVector Boost(const ThreeVector&) noexcept;

        /// Boost the four vector to the frame given by the three vector velocities
        ///@param beta_x: The boost in the x-direction
        ///@param beta_y: The boost in the y-direction
        ///@param beta_z: The boost in the z-direction
        ///@return FourVector: The vector in the corresponding frame
        FourVector Boost(const double&, const double&, const double&) noexcept;

        /// Calculate the cross product between two four vectors
        ///@param other: The vector to take the cross product with respect to
        ///@return FourVector: The vector perpendicular to the two inputs
        FourVector Cross(const FourVector&) const noexcept;

        /// Return the the boost vector required to boost to the rest frame of a given
        /// four vector
        ///@return ThreeVector: Boost vector
        ThreeVector BoostVector() const noexcept;

        /// Calculate the dot product between two four vectors
        ///@param other: four vector to take dot product with
        ///@return double: the value of the dot product
        double Dot(const FourVector& other) const noexcept {
            return (*this) * other;
        }

        /// Return a string representation of the vector
        ///@return std::string: a string representation of the vector
        std::string ToString() const noexcept;
        ///@}

        /// @name Operator Overloads
        /// @{
        /// Operator overloads of math functions for ease of use

        /// Assign a vector to another vector object
        ///@param other: four vector to be assigned
        ///@return FourVector: assigned new vector
        FourVector& operator=(const FourVector&) noexcept = default;
        FourVector& operator=(FourVector&&) noexcept = default;

        /// Add two four vectors together
        ///@param other: four vector to add to this one
        ///@return FourVector: The sum of the two four vectors
        FourVector& operator+=(const FourVector&) noexcept;

        /// Subtract two four vectors together
        ///@param other: four vector to subtract from this one
        ///@return FourVector: The difference of the two four vectors
        FourVector& operator-=(const FourVector&) noexcept;

        /// Scale a four vector by a constant
        ///@param scale: value to scale four vector by
        ///@return FourVector: The scaled four vector
        FourVector& operator*=(const double&) noexcept;

        /// Reduce the magnitude of a four vector by a constant
        ///@param scale: value to scale four vector by
        ///@return FourVector: The scaled four vector
        FourVector& operator/=(const double&);

        /// Calculate dot product of two vectors
        ///@param other: The other vector to take the dot product with
        ///@return double: The value of the dot product
        double operator*(const FourVector&) const noexcept;

        /// Negate a given vector
        ///@return FourVector: The negative of the input vector
        FourVector operator-() const noexcept;

        /// Unary plus operator
        ///@return FourVector: The input vector
        FourVector operator+() const noexcept;

        /// Scale a four vector by a constant
        ///@param scale: value to scale four vector by
        ///@return FourVector: The scaled four vector
        FourVector operator*(const double&) const noexcept;

        /// Reduce the magnitude of a four vector by a constant
        ///@param scale: value to scale four vector by
        ///@return FourVector: The scaled four vector
        FourVector operator/(const double&) const;

        /// Calculate difference of two vectors
        ///@param other: The other vector to take the difference with
        ///@return FourVector: The resulting difference vector
        FourVector operator-(const FourVector&) const noexcept;

        /// Calculate sum of two vectors
        ///@param other: The other vector to take the sum with
        ///@return FourVector: The resulting sum vector
        FourVector operator+(const FourVector&) const noexcept;

        // Comparison Operators

        /// Determine if two four vectors are equivalent
        ///@param other: The vector to compare against
        ///@return bool: True, if the vectors are equal otherwise False
        bool operator==(const FourVector&) const noexcept;

        /// Determine if two four vectors are not equivalent
        ///@param other: The vector to compare against
        ///@return bool: False, if the vectors are equal otherwise True
        bool operator!=(const FourVector& other) const noexcept {return !(*this == other);}

        // Access Operators

        /// Access a given index from the vector
        ///@param idx: Index to access
        ///@return double: The value of the vector at the given index
        double& operator[](const std::size_t& idx) {
            if(idx > 3) throw std::range_error("Max value is 3.");
            return vec[idx];
        }

        /// Access a given index from the vector
        ///@param idx: Index to access
        ///@return double: The value of the vector at the given index
        const double& operator[](const std::size_t& idx) const {
            if(idx > 3) throw std::range_error("Max value is 3.");
            return vec[idx];
        }
        /// @}

        ///@name Stream Operators
        /// @{
        /// Stream operators for writing to and reading from streams

        /// Write out a vector to an output stream
        ///@param ostr: Output stream to write to
        ///@param vec: The four vector to be written out
        template<typename OStream>
        friend OStream& operator<<(OStream &os, const FourVector &vec4) {
            os << "FourVector(" << vec4.Px() << ", " << vec4.Py() << ", " << vec4.Pz() << ", " << vec4.E() << ")";
            return os;
        }

        /// Write in a vector to an input stream
        ///@param istr: Input stream to read from
        ///@param vec: The four vector to be read into
        friend std::istream& operator>>(std::istream&, FourVector&);
        /// @}

    private:
        std::array<double, 4> vec;
};

FourVector operator*(const double&, const FourVector&) noexcept;

}

#endif /* end of include guard: FOURVECTOR_HH */
