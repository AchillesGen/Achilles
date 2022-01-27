#ifndef PLUGINS_CHANNELS1_HH
#define PLUGINS_CHANNELS1_HH

// Channel_Generator_UniV
#include "PHASIC++/Channels/Channel_Elements.H"
#include "plugins/Sherpa/Channels.hh"
#include "nuchic/PhaseSpaceFactory.hh"

using namespace PHASIC;
using namespace ATOOLS;

namespace PHASIC {

class Channels1 : public Channels {
    protected:
        constexpr static double m_amct=1.,m_alpha=0.7,m_ctmax=1.0,m_ctmin=-1.;
        double s2, s3;

    private:
        std::string m_name;

    public:
        Channels1(const std::vector<double> &s, std::string name) 
              : s2{s[0]}, s3{s[1]}, m_name{std::move(name)} {}

        double GenerateWeight(const std::vector<Vec4D>&, std::vector<double>&) override = 0;
        void GeneratePoint(std::vector<Vec4D>&, const std::vector<double>&) override = 0;
        size_t NDims() const override { return 2; }
        YAML::Node ToYAML() const override {
            YAML::Node result;
            result["Name"] = m_name;
            result["Masses"] = std::vector<double>{s2, s3};
            return result;
        }
};

class C1_0 : public Channels1, nuchic::RegistrablePS<Channels, C1_0, std::vector<double>> {
    public:
        C1_0(const std::vector<double> &s) : Channels1(s, Name()) {}
        static std::string Name() { return "C1_0"; }
        static std::unique_ptr<Channels> Construct(const std::vector<double> &s) {
            if(s.size() != 2) {
                auto msg = fmt::format("Incorrect number of masses. Expected 2. Got {}", s.size());
                throw std::runtime_error(msg);
            }
            return std::make_unique<C1_0>(s);
        }

        double GenerateWeight(const std::vector<Vec4D>&, std::vector<double> &) override;
        void GeneratePoint(std::vector<Vec4D>&, const std::vector<double> &) override;
};

class C1_1 : public Channels1, nuchic::RegistrablePS<Channels, C1_1, std::vector<double>> {
    public:
        C1_1(const std::vector<double> &s) : Channels1(s, Name()) {}
        static std::string Name() { return "C1_1"; }
        static std::unique_ptr<Channels> Construct(const std::vector<double> &s) {
            if(s.size() != 2) {
                auto msg = fmt::format("Incorrect number of masses. Expected 2. Got {}", s.size());
                throw std::runtime_error(msg);
            }
            return std::make_unique<C1_1>(s);
        }
        double GenerateWeight(const std::vector<Vec4D>&, std::vector<double> &) override;
        void GeneratePoint(std::vector<Vec4D>&, const std::vector<double> &) override;
};

class C1_2 : public Channels1, nuchic::RegistrablePS<Channels, C1_2, std::vector<double>> {
    public:
        C1_2(const std::vector<double> &s) : Channels1(s, Name()) {}
        static std::string Name() { return "C1_2"; }
        static std::unique_ptr<Channels> Construct(const std::vector<double> &s) {
            if(s.size() != 2) {
                auto msg = fmt::format("Incorrect number of masses. Expected 2. Got {}", s.size());
                throw std::runtime_error(msg);
            }
            return std::make_unique<C1_2>(s);
        }

        double GenerateWeight(const std::vector<Vec4D>&, std::vector<double> &) override;
        void GeneratePoint(std::vector<Vec4D>&, const std::vector<double> &) override;
};

}

#endif
