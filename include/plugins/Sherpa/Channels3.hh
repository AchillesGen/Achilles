#ifndef PLUGINS_CHANNELS3_HH
#define PLUGINS_CHANNELS3_HH

// Channel_Generator_UniV
#include "PHASIC++/Channels/Channel_Elements.H"
#include "plugins/Sherpa/Channels.hh"

using namespace PHASIC;
using namespace ATOOLS;

namespace PHASIC {

class Channels3 : public Channels {
    protected:
        constexpr static double m_salpha=0.9,m_thexp=0.9;
        constexpr static double m_amct=1.,m_alpha=0.9,m_ctmax=1.,m_ctmin=-1.;
        double s2, s3, s4, s5;

        void GenerateNuclearPoint(std::vector<Vec4D>&, const std::vector<double>&) const;
        double GenerateNuclearWeight(const std::vector<Vec4D>&, std::vector<double>&) const;

    private:
        std::string m_name;

    public:
        Channels3(const std::vector<double> &s, std::string name) 
            : s2{std::move(s[0])}, s3{std::move(s[1])}, s4{std::move(s[2])}, s5{std::move(s[3])},
              m_name{std::move(name)} {}

        double GenerateWeight(const std::vector<Vec4D>&, std::vector<double>&) override = 0;
        void GeneratePoint(std::vector<Vec4D>&, const std::vector<double>&) override = 0;
        size_t NDims() const override { return 5; }
        YAML::Node ToYAML() const override {
            YAML::Node result;
            result["Name"] = m_name;
            result["Masses"] = std::vector<double>{s2, s3, s4, s5};
            return result;
        }
};

class C3_0 : public Channels3, achilles::RegistrablePS<Channels, C3_0, std::vector<double>> {
    public:
        C3_0(const std::vector<double> &s) : Channels3(s, Name()) {}
        static std::string Name() { return "C3_0"; }
        static std::unique_ptr<Channels> Construct(const std::vector<double> &s) {
            if(s.size() != 4) {
                auto msg = fmt::format("Incorrect number of masses. Expected 4. Got {}", s.size());
                throw std::runtime_error(msg);
            }
            return std::make_unique<C3_0>(s);
        }

        double GenerateWeight(const std::vector<Vec4D>&, std::vector<double> &) override;
        void GeneratePoint(std::vector<Vec4D>&, const std::vector<double> &) override;
};

class C3_1 : public Channels3, achilles::RegistrablePS<Channels, C3_1, std::vector<double>> {
    public:
        C3_1(const std::vector<double> &s) : Channels3(s, Name()) {}
        static std::string Name() { return "C3_1"; }
        static std::unique_ptr<Channels> Construct(const std::vector<double> &s) {
            if(s.size() != 4) {
                auto msg = fmt::format("Incorrect number of masses. Expected 4. Got {}", s.size());
                throw std::runtime_error(msg);
            }
            return std::make_unique<C3_1>(s);
        }

        double GenerateWeight(const std::vector<Vec4D>&, std::vector<double> &) override;
        void GeneratePoint(std::vector<Vec4D>&, const std::vector<double> &) override;
};

class C3_2 : public Channels3, achilles::RegistrablePS<Channels, C3_2, std::vector<double>> {
    public:
        C3_2(const std::vector<double> &s) : Channels3(s, Name()) {}
        static std::string Name() { return "C3_2"; }
        static std::unique_ptr<Channels> Construct(const std::vector<double> &s) {
            if(s.size() != 4) {
                auto msg = fmt::format("Incorrect number of masses. Expected 4. Got {}", s.size());
                throw std::runtime_error(msg);
            }
            return std::make_unique<C3_2>(s);
        }

        double GenerateWeight(const std::vector<Vec4D>&, std::vector<double> &) override;
        void GeneratePoint(std::vector<Vec4D>&, const std::vector<double> &) override;
};

class C3_3 : public Channels3, achilles::RegistrablePS<Channels, C3_3, std::vector<double>> {
    public:
        C3_3(const std::vector<double> &s) : Channels3(s, Name()) {}
        static std::string Name() { return "C3_3"; }
        static std::unique_ptr<Channels> Construct(const std::vector<double> &s) {
            if(s.size() != 4) {
                auto msg = fmt::format("Incorrect number of masses. Expected 4. Got {}", s.size());
                throw std::runtime_error(msg);
            }
            return std::make_unique<C3_3>(s);
        }

        double GenerateWeight(const std::vector<Vec4D>&, std::vector<double> &) override;
        void GeneratePoint(std::vector<Vec4D>&, const std::vector<double> &) override;
};

class C3_4 : public Channels3, achilles::RegistrablePS<Channels, C3_4, std::vector<double>> {
    public:
        C3_4(const std::vector<double> &s) : Channels3(s, Name()) {}
        static std::string Name() { return "C3_4"; }
        static std::unique_ptr<Channels> Construct(const std::vector<double> &s) {
            if(s.size() != 4) {
                auto msg = fmt::format("Incorrect number of masses. Expected 4. Got {}", s.size());
                throw std::runtime_error(msg);
            }
            return std::make_unique<C3_4>(s);
        }

        double GenerateWeight(const std::vector<Vec4D>&, std::vector<double> &) override;
        void GeneratePoint(std::vector<Vec4D>&, const std::vector<double> &) override;
};

class C3_5 : public Channels3, achilles::RegistrablePS<Channels, C3_5, std::vector<double>> {
    public:
        C3_5(const std::vector<double> &s) : Channels3(s, Name()) {}
        static std::string Name() { return "C3_5"; }
        static std::unique_ptr<Channels> Construct(const std::vector<double> &s) {
            if(s.size() != 4) {
                auto msg = fmt::format("Incorrect number of masses. Expected 4. Got {}", s.size());
                throw std::runtime_error(msg);
            }
            return std::make_unique<C3_5>(s);
        }

        double GenerateWeight(const std::vector<Vec4D>&, std::vector<double> &) override;
        void GeneratePoint(std::vector<Vec4D>&, const std::vector<double> &) override;
};

class C3_6 : public Channels3, achilles::RegistrablePS<Channels, C3_6, std::vector<double>> {
    public:
        C3_6(const std::vector<double> &s) : Channels3(s, Name()) {}
        static std::string Name() { return "C3_6"; }
        static std::unique_ptr<Channels> Construct(const std::vector<double> &s) {
            if(s.size() != 4) {
                auto msg = fmt::format("Incorrect number of masses. Expected 4. Got {}", s.size());
                throw std::runtime_error(msg);
            }
            return std::make_unique<C3_6>(s);
        }

        double GenerateWeight(const std::vector<Vec4D>&, std::vector<double> &) override;
        void GeneratePoint(std::vector<Vec4D>&, const std::vector<double> &) override;
};

class C3_7 : public Channels3, achilles::RegistrablePS<Channels, C3_7, std::vector<double>> {
    public:
        C3_7(const std::vector<double> &s) : Channels3(s, Name()) {}
        static std::string Name() { return "C3_7"; }
        static std::unique_ptr<Channels> Construct(const std::vector<double> &s) {
            if(s.size() != 4) {
                auto msg = fmt::format("Incorrect number of masses. Expected 4. Got {}", s.size());
                throw std::runtime_error(msg);
            }
            return std::make_unique<C3_7>(s);
        }

        double GenerateWeight(const std::vector<Vec4D>&, std::vector<double> &) override;
        void GeneratePoint(std::vector<Vec4D>&, const std::vector<double> &) override;
};

}

#endif
