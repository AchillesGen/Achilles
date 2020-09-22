#include <iostream>

#include "nuchic/FourVector.hh"
#include "nuchic/ThreeVector.hh"

using namespace nuchic;

FourVector::FourVector(const ThreeVector& other, const double& E) noexcept
    : vec({other[0], other[1], other[2], E}) {}

double FourVector::M() const noexcept {
    return std::sqrt(M2());
}

double FourVector::Theta() const noexcept {
    return atan2(Pt(), Pz());
}

double FourVector::Phi() const noexcept {
    const double phi = atan2(Py(), Px());
    if(phi < 0) return phi + 2*M_PI;
    return phi;
}

double FourVector::Rapidity() const noexcept {
    return log((E() + Pz())/(E() - Pz()))/2;
}

double FourVector::DeltaR(const FourVector& other) const noexcept {
    const double DEta = Rapidity() - other.Rapidity();
    const double DPhi = Phi() - other.Phi();
    return sqrt(DEta*DEta + DPhi*DPhi);
}

ThreeVector FourVector::Vec3() const noexcept {
    return {vec[0], vec[1], vec[2]};
}

void FourVector::SetVectM(const ThreeVector& vec3, const double& mass) noexcept {
    vec[0] = vec3.Px();
    vec[1] = vec3.Py();
    vec[2] = vec3.Pz();
    vec[3] = sqrt(mass*mass + vec3*vec3);
}

FourVector FourVector::Boost(const ThreeVector& beta) noexcept {
    const double beta2 = beta*beta;
    const double gamma = 1.0/sqrt(1.0 - beta2);
    const double betap = beta[0]*vec[0] + beta[1]*vec[1] + beta[2]*vec[2];
    const double gamma2 = beta2 > 0 ? (gamma-1.0)/beta2 : 0.0;

    double pX = vec[0] + gamma2*betap*beta[0]+gamma*beta[0]*vec[3];
    double pY = vec[1] + gamma2*betap*beta[1]+gamma*beta[1]*vec[3];
    double pZ = vec[2] + gamma2*betap*beta[2]+gamma*beta[2]*vec[3];
    double E = gamma*(vec[3] + betap);

    return {pX, pY, pZ, E};
}

FourVector FourVector::Boost(const double& beta_x, const double& beta_y,
                                             const double& beta_z) noexcept {
    return Boost(ThreeVector(beta_x, beta_y, beta_z));
}

FourVector FourVector::Cross(const FourVector& other) const noexcept {
    ThreeVector result = this -> Vec3().Cross(other.Vec3());
    return {result, 0};
}

ThreeVector FourVector::BoostVector() const noexcept {
    return {vec[0]/vec[3], vec[1]/ vec[3], vec[2]/vec[3]};
}

FourVector& FourVector::operator+=(const FourVector& other) noexcept {
    vec[0] += other.vec[0];
    vec[1] += other.vec[1];
    vec[2] += other.vec[2];
    vec[3] += other.vec[3];

    return *this;
}

FourVector& FourVector::operator-=(const FourVector& other) noexcept {
    vec[0] -= other.vec[0];
    vec[1] -= other.vec[1];
    vec[2] -= other.vec[2];
    vec[3] -= other.vec[3];

    return *this;
}

FourVector& FourVector::operator*=(const double& scale) noexcept {
    vec[0] *= scale;
    vec[1] *= scale;
    vec[2] *= scale;
    vec[3] *= scale;

    return *this;
}

FourVector& FourVector::operator/=(const double& scale) {
    vec[0] /= scale;
    vec[1] /= scale;
    vec[2] /= scale;
    vec[3] /= scale;

    return *this;
}

double FourVector::operator*(const FourVector& other) const noexcept {
    return vec[3]*other.vec[3] 
        - (vec[0]*other.vec[0] + vec[1]*other.vec[1] + vec[2]*other.vec[2]);
}

FourVector FourVector::operator-() const noexcept {
   return {-vec[0],-vec[1],-vec[2],-vec[3]};
}

FourVector FourVector::operator+() const noexcept {
    return *this;
}

FourVector FourVector::operator*(const double& scale) const noexcept {
    return FourVector(*this)*=scale;
}

FourVector FourVector::operator/(const double& scale) const {
    return FourVector(*this)/=scale;
}

FourVector FourVector::operator+(const FourVector& other) const noexcept {
    return FourVector(*this)+=other;
}

FourVector FourVector::operator-(const FourVector& other) const noexcept {
    return FourVector(*this)-=other;
}

bool FourVector::operator==(const FourVector& other) const noexcept {
    return vec==other.vec;
}

std::string FourVector::ToString() const noexcept {
    return "FourVector(" + std::to_string(vec[0]) + ", " + std::to_string(vec[1])
        + ", " + std::to_string(vec[2]) + ", " + std::to_string(vec[3]) + ")";
}

namespace nuchic {

std::istream& operator>>(std::istream& is, FourVector& vec) {
    std::string head(11, ' '), sep1(1, ' '), sep2(1, ' '),
        sep3(1, ' '), tail(1, ' ');
    double px, py, pz, e;
    is.read(&head[0], 11);
    is >> px;
    is.read(&sep1[0], 1);
    is >> py;
    is.read(&sep2[0], 1);
    is >> pz;
    is.read(&sep3[0], 1);
    is >> e;
    is.read(&tail[0], 1);
    if(head == "FourVector(" &&
       sep1 == "," && sep2 == "," && 
       sep3 == "," && tail == ")") 
        vec = FourVector(px, py, pz, e);
    return is;
}

FourVector operator*(const double& s, const FourVector& v) noexcept {
    return v*s;
}

}
