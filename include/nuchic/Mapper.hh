#ifndef MAPPER_HH
#define MAPPER_HH

#include <memory>
#include <string>
#include <vector>

namespace nuchic {

template<typename T>
class Mapper {
    public:
        template<typename C>
        using Mapper_ptr = std::unique_ptr<Mapper<C>>;
        template<typename C>
        using Mapper_sptr = std::shared_ptr<Mapper<C>>;

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
        virtual size_t NDims() const = 0;

    private:
        std::string mapping_name{};
};


}

#endif
