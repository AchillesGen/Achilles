#ifndef PHASESPACE_MAPPER_HH
#define PHASESPACE_MAPPER_HH

#include "nuchic/Mapper.hh"

namespace nuchic {

class FourVector;

class PSMapper : public Mapper<FourVector> {
    public:
        PSMapper(size_t _nleptons, size_t _nhadrons,
                 Mapper_sptr<FourVector> _initial, Mapper_ptr<FourVector> _leptonic,
                 Mapper_sptr<FourVector> _hadronic)
            : nleptons{std::move(_nleptons)}, nhadrons{std::move(_nhadrons)},
              initial{std::move(_initial)}, leptonic{std::move(_leptonic)},
              hadronic{std::move(_hadronic)} {}

        void GeneratePoint(std::vector<FourVector>&, const std::vector<double>&) const override;
        double GenerateWeight(const std::vector<FourVector>&, std::vector<double>&) const override;
        size_t NDims() const override { 
            return initial -> NDims() + leptonic -> NDims() + hadronic -> NDims();
        }

    private:
        size_t nleptons, nhadrons;
        Mapper_sptr<FourVector> initial;
        Mapper_ptr<FourVector> leptonic; 
        Mapper_sptr<FourVector> hadronic;
};

}

#endif
