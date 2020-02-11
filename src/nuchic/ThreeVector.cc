#include <iostream>

#include "nuchic/ThreeVector.hh"

const double nuchic::ThreeVector::Theta() const noexcept {
    return atan2(Pt(), vec[3]);
}

const double nuchic::ThreeVector::Phi() const noexcept {
    const double phi = atan2(vec[1], vec[0]);
    if(phi < 0) return phi + 2*M_PI;
    return phi;
}

const nuchic::ThreeVector nuchic::ThreeVector::Cross(const nuchic::ThreeVector& other) const noexcept {
    const double px = vec[1]*other.vec[2] - vec[2]*other.vec[1];    
    const double py = vec[2]*other.vec[0] - vec[0]*other.vec[2];
    const double pz = vec[0]*other.vec[1] - vec[1]*other.vec[0];    

    return nuchic::ThreeVector(px, py, pz);
}

const nuchic::ThreeVector nuchic::ThreeVector::Unit() const noexcept {
    const double norm = Magnitude();

    return nuchic::ThreeVector(vec[0]/norm, vec[1]/norm, vec[2]/norm);
}

nuchic::ThreeVector& nuchic::ThreeVector::operator+=(const nuchic::ThreeVector& other) noexcept {
    vec[0] += other.vec[0];
    vec[1] += other.vec[1];
    vec[2] += other.vec[2];

    return *this;
}

nuchic::ThreeVector& nuchic::ThreeVector::operator-=(const nuchic::ThreeVector& other) noexcept {
    vec[0] -= other.vec[0];
    vec[1] -= other.vec[1];
    vec[2] -= other.vec[2];

    return *this;
}

nuchic::ThreeVector& nuchic::ThreeVector::operator*=(const double& scale) noexcept {
    vec[0] *= scale;
    vec[1] *= scale;
    vec[2] *= scale;

    return *this;
}

nuchic::ThreeVector& nuchic::ThreeVector::operator/=(const double& scale) {
    vec[0] /= scale;
    vec[1] /= scale;
    vec[2] /= scale;

    return *this;
}

double nuchic::ThreeVector::operator*(const nuchic::ThreeVector& other) const noexcept {
    return vec[0]*other.vec[0] + vec[1]*other.vec[1] + vec[2]*other.vec[2];
}

nuchic::ThreeVector nuchic::ThreeVector::operator-() const noexcept {
   return nuchic::ThreeVector(-vec[0], -vec[1], -vec[2]);
}

nuchic::ThreeVector nuchic::ThreeVector::operator+() const noexcept {
    return nuchic::ThreeVector(*this);
}

nuchic::ThreeVector nuchic::ThreeVector::operator*(const double& scale) const noexcept {
    return nuchic::ThreeVector(*this)*=scale;
}

nuchic::ThreeVector nuchic::ThreeVector::operator/(const double& scale) const {
    return nuchic::ThreeVector(*this)/=scale;
}

nuchic::ThreeVector nuchic::ThreeVector::operator+(const ThreeVector& other) const noexcept {
    return nuchic::ThreeVector(*this)+=other;
}

nuchic::ThreeVector nuchic::ThreeVector::operator-(const ThreeVector& other) const noexcept {
    return nuchic::ThreeVector(*this)-=other;
}

bool nuchic::ThreeVector::operator==(const nuchic::ThreeVector& other) const noexcept {
    return vec==other.vec;
}

const std::string nuchic::ThreeVector::ToString() const noexcept {
    return "ThreeVector(" + std::to_string(vec[0]) + ", " 
        + std::to_string(vec[1]) + ", " + std::to_string(vec[2]) + ")";
}

std::ostream& operator<<(std::ostream& os, const nuchic::ThreeVector& vec) {
    os << "ThreeVector(" << vec.Px() << ", " << vec.Py() << ", " << vec.Pz() << ")";
    return os;
}

std::istream& operator>>(std::istream& is, nuchic::ThreeVector& vec) {
    std::string head(12, ' '), sep1(1, ' '), sep2(1, ' '), tail(1, ' ');
    double px, py, pz, e;
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
        vec = nuchic::ThreeVector(px, py, pz);
    return is;
}

nuchic::ThreeVector operator*(const double& s, const nuchic::ThreeVector& v) noexcept {
    return v*s;
}
