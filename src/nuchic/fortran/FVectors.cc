#include "nuchic/FourVector.hh"
#include "nuchic/ThreeVector.hh"

extern "C" {

    using Vec4 = nuchic::FourVector;
    using Vec3 = nuchic::ThreeVector;

    // FourVector Constructor and Destructor
    Vec4* CreateFourVector(const double x, const double y, const double z, const double E) {
        return new Vec4(x, y, z, E);
    }
    void DeleteFourVector(Vec4 *self) {
        delete self;
    }

    // Functions with FourVectors
    Vec4* Boost(const Vec4* self, const Vec3* beta) {
        return new Vec4(self -> Boost(*beta));
    }
    Vec3* BoostVector(const Vec4* self) {
        return new Vec3(self -> BoostVector());
    }
    double Dot4(const Vec4* self, const Vec4* other) {
        return self -> Dot(*other);
    }
    Vec4* Add4(const Vec4* self, const Vec4* other) {
        return new Vec4((*self) + (*other));
    }
    Vec4* Sub4(const Vec4* self, const Vec4* other) {
        return new Vec4((*self) - (*other));
    }
    Vec4* Scale4(const Vec4* self, const double scale) {
        return new Vec4(scale*(*self));
    }
    double Get4(const Vec4* self, const int idx) {
        return self -> operator[](static_cast<size_t>(idx));
    }
    void Print4(const Vec4* self) {
        fmt::print("{}\n", *self);
    }

    // ThreeVector Constructor and Destructor
    Vec3* CreateThreeVector(const double x, const double y, const double z) {
        return new Vec3(x, y, z);
    }
    Vec3* New3(Vec3 *self) {
        return new Vec3(*self);
    }

    void DeleteThreeVector(Vec3 *self) {
        delete self;
    }

    // Functions with ThreeVectors
    double Dot3(const Vec3* self, const Vec3* other) {
        return self -> Dot(*other);
    }
    Vec3* Add3(const Vec3* self, const Vec3* other) {
        return new Vec3((*self) + (*other));
    }
    Vec3* Sub3(const Vec3* self, const Vec3* other) {
        return new Vec3((*self) - (*other));
    }
    Vec3* Scale3(const Vec3* self, const double scale) {
        return new Vec3(scale*(*self));
    }
    double Get3(const Vec3* self, const int idx) {
        return self -> operator[](static_cast<size_t>(idx));
    }
    void Print3(const Vec3* self) {
        fmt::print("{}\n", *self);
    }

}
