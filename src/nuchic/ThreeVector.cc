#include <iostream>

#include "nuchic/ThreeVector.hh"

using namespace nuchic;

double ThreeVector::Theta() const noexcept {
    return atan2(Pt(), Pz());
}

double ThreeVector::Phi() const noexcept {
    const double phi = atan2(Py(), Px());
    if(phi < 0) return phi + 2*M_PI;
    return phi;
}

ThreeVector ThreeVector::Cross(const ThreeVector& other) const noexcept {
    const double px = vec[1]*other.vec[2] - vec[2]*other.vec[1];    
    const double py = vec[2]*other.vec[0] - vec[0]*other.vec[2];
    const double pz = vec[0]*other.vec[1] - vec[1]*other.vec[0];    

    return {px, py, pz};
}

ThreeVector ThreeVector::Unit() const noexcept {
    const double norm = Magnitude();

    return {vec[0]/norm, vec[1]/norm, vec[2]/norm};
}

ThreeVector& ThreeVector::operator+=(const ThreeVector& other) noexcept {
    vec[0] += other.vec[0];
    vec[1] += other.vec[1];
    vec[2] += other.vec[2];

    return *this;
}

ThreeVector& ThreeVector::operator-=(const ThreeVector& other) noexcept {
    vec[0] -= other.vec[0];
    vec[1] -= other.vec[1];
    vec[2] -= other.vec[2];

    return *this;
}

ThreeVector& ThreeVector::operator*=(const double& scale) noexcept {
    vec[0] *= scale;
    vec[1] *= scale;
    vec[2] *= scale;

    return *this;
}

ThreeVector& ThreeVector::operator/=(const double& scale) {
    vec[0] /= scale;
    vec[1] /= scale;
    vec[2] /= scale;

    return *this;
}

double ThreeVector::operator*(const ThreeVector& other) const noexcept {
    return vec[0]*other.vec[0] + vec[1]*other.vec[1] + vec[2]*other.vec[2];
}

ThreeVector ThreeVector::operator-() const noexcept {
   return {-vec[0], -vec[1], -vec[2]};
}

ThreeVector ThreeVector::operator+() const noexcept {
    return *this;
}

ThreeVector ThreeVector::operator*(const double& scale) const noexcept {
    return ThreeVector(*this)*=scale;
}

ThreeVector ThreeVector::operator/(const double& scale) const {
    return ThreeVector(*this)/=scale;
}

ThreeVector ThreeVector::operator+(const ThreeVector& other) const noexcept {
    return ThreeVector(*this)+=other;
}

ThreeVector ThreeVector::operator-(const ThreeVector& other) const noexcept {
    return ThreeVector(*this)-=other;
}

bool ThreeVector::operator==(const ThreeVector& other) const noexcept {
    return vec==other.vec;
}

std::string ThreeVector::ToString() const noexcept {
    return "ThreeVector(" + std::to_string(vec[0]) + ", " 
        + std::to_string(vec[1]) + ", " + std::to_string(vec[2]) + ")";
}

namespace nuchic {

std::istream& operator>>(std::istream& is, ThreeVector& vec) {
    std::string head(12, ' '), sep1(1, ' '), sep2(1, ' '), tail(1, ' ');
    double px, py, pz;
    is.read(&head[0], 12);
    is >> px;
    is.read(&sep1[0], 1);
    is >> py;
    is.read(&sep2[0], 1);
    is >> pz;
    is.read(&tail[0], 1);
    if(head == "ThreeVector(" &&
       sep1 == "," && sep2 == "," && 
       tail == ")") 
        vec = ThreeVector(px, py, pz);
    return is;
}

ThreeVector operator*(const double& s, const ThreeVector& v) noexcept {
    return v*s;
}

}
