#ifndef MAPPER_HH
#define MAPPER_HH

#include <string>
#include <vector>

namespace nuchic {

template<typename T>
class Mapper {
    public:
        Mapper() = default;
        Mapper(std::string name) : mapping_name{std::move(name)} {}
        Mapper(const Mapper&) = delete;
        Mapper(Mapper&&) = delete;
        Mapper& operator=(const Mapper&) = delete;
        Mapper& operator=(Mapper&&) = delete;
        virtual ~Mapper() = default;

        // Functions
        virtual void GeneratePoint(std::vector<T>&, const std::vector<double> &) const = 0;
        virtual double GenerateWeight(const std::vector<T>&, std::vector<double> &) const = 0;

    private:
        std::string mapping_name{};
};

}

#endif
