#ifndef RAMBO_HH
#define RAMBO_HH

#include "nuchic/FourVector.hh"

namespace nuchic {

class Rambo {
    private:
        static constexpr size_t nin = 2;
        static constexpr size_t NMAX = 20;

    public:
        Rambo(size_t);

        void GeneratePoint(std::vector<FourVector>&,
                           const std::vector<double>&,
                           const std::vector<double>&);
        void GenerateWeight(const std::vector<FourVector>&,
                            const std::vector<double>&);
        double Weight() const { return weight; }
        size_t Nin() const { return nin; }
        size_t Nout() const { return nout; }
        void Nout(size_t nout_) { nout = nout_; }
        size_t Npart() const { return nin + nout; }

    private:
        void MassivePoint(std::vector<FourVector>&, 
                          const std::vector<double>&,
                          double) const;
        void MassiveWeight(const std::vector<FourVector>&,
                           const std::vector<double>&,
                           double);

        std::array<double, NMAX-nin+1> Z{};

        size_t nout;
        double weight{};
};

}

#endif
