#ifndef PARTICLE_HH
#define PARTICLE_HH

#include <iosfwd>
#include <vector>

#include "nuchic/ThreeVector.hh"
#include "nuchic/FourVector.hh"

class Particle {
    public:
        // Constructors and Destructors
        Particle(const int& _pid = 0, const FourVector& mom = FourVector(),
                 const ThreeVector& pos = ThreeVector(), const int& _status = 0,
                 const std::vector<int>& _mothers = std::vector<int>(),
                 const std::vector<int>& _daughters = std::vector<int>()) noexcept :
            pid(_pid), momentum(mom), position(pos), status(_status),
            mothers(_mothers), daughters(_daughters) {formationZone = 0;}
        ~Particle() {};

        // Setters
        void SetPID(const int& _pid) noexcept {pid = _pid;}
        void SetPosition(const ThreeVector& pos) noexcept {position = pos;}
        void SetMomentum(const FourVector& mom) noexcept {momentum = mom;}
        void SetStatus(const int& _status) noexcept {status = _status;}
        void SetMothers(const std::vector<int>& _mothers) noexcept {mothers = _mothers;}
        void SetDaughters(const std::vector<int>& _daughters) noexcept {daughters = _daughters;}
        void AddMother(const int& idx) noexcept {mothers.push_back(idx);}
        void AddDaughter(const int& idx) noexcept {daughters.push_back(idx);}
        void SetFormationZone(const FourVector&, const FourVector&) noexcept;
        void UpdateFormationZone(const double& timeStep) noexcept {formationZone -= timeStep;}

        // Getters
        const int& PID() const noexcept {return pid;}
        const ThreeVector& Position() const noexcept {return position;}
        const FourVector& Momentum() const noexcept {return momentum;}
        const ThreeVector Beta() const noexcept {return momentum.BoostVector();}
        const int& Status() const noexcept {return status;}
        const std::vector<int>& Mothers() const noexcept {return mothers;}
        const std::vector<int>& Daughters() const noexcept {return daughters;}
        const double& FormationZone() const noexcept {return formationZone;}
        const double Mass() const noexcept {return momentum.M();}
        const double Px() const noexcept {return momentum.Px();}
        const double Py() const noexcept {return momentum.Py();}
        const double Pz() const noexcept {return momentum.Pz();}
        const double E() const noexcept {return momentum.E();}
        const double Radius() const noexcept {return position.Magnitude();}

        // Logical Functions
        const bool InFormationZone() const noexcept {return formationZone > 0;}
        const bool IsBackground() const noexcept {return status == 0;}
        const bool IsPropagating() const noexcept {return status == -1;}
        const bool IsFinal() const noexcept {return status == 1;}

        // Functions
        void Propagate(const double&) noexcept;
        void BackPropagate(const double&) noexcept;
        const std::string ToString() const noexcept;

        // Comparison Operators
        bool operator==(const Particle&) const noexcept;
        bool operator!=(const Particle& other) const noexcept {return !(*this == other);}
        
        // Stream Operators
        friend std::ostream& operator<<(std::ostream&, const Particle&);
        friend std::istream& operator>>(std::istream&, Particle&);

    private:
        int pid, status;
        std::vector<int> mothers, daughters;
        double formationZone;
        ThreeVector position;
        FourVector momentum;
};

#endif // end of include guard: PARTICLE_HH
