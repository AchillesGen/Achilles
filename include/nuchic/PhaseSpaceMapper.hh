#ifndef PHASESPACE_MAPPER_HH
#define PHASESPACE_MAPPER_HH

#include "nuchic/Mapper.hh"

namespace nuchic {

class FourVector;

class PSMapper : public Mapper<FourVector> {
    public:
        PSMapper(size_t _nleptons, size_t _nhadrons,
                 Mapper_sptr<FourVector> _lbeam, Mapper_sptr<FourVector> _hbeam,
                 Mapper_ptr<FourVector> _main)
            : nleptons{std::move(_nleptons)}, nhadrons{std::move(_nhadrons)},
              lbeam{std::move(_lbeam)},
              hbeam{std::move(_hbeam)},
              main{std::move(_main)} {}

        void GeneratePoint(std::vector<FourVector>&, const std::vector<double>&) const override;
        double GenerateWeight(const std::vector<FourVector>&, std::vector<double>&) const override;
        size_t NDims() const override { 
            return lbeam -> NDims() + hbeam -> NDims() + main -> NDims();
        }

    private:
        size_t nleptons, nhadrons;
        Mapper_sptr<FourVector> lbeam, hbeam;
        Mapper_ptr<FourVector> main; 
};

}

#endif
