#ifndef INTERACTIONS_HH
#define INTERACTIONS_HH

#include <array>
#include <map>
#include <vector>

#include "H5Cpp.h"

#include "nuchic/Interpolation.hh"

class ThreeVector;
class Particle;

double CrossSection(bool, const double&);
double CrossSectionAngle(bool, const double&, const double&);
double CrossSectionLab(bool, const double&) noexcept;
ThreeVector MakeMomentumAngular(bool, const double&, const double&, const std::array<double, 2>&);

class Interactions {
    public:
        // Constructors
        Interactions() {}
        virtual ~Interactions() {}

        // Functions
        virtual bool IsRegistered() const noexcept = 0;
        virtual double CrossSection(const Particle&, const Particle&) const = 0;
        virtual ThreeVector MakeMomentum(bool, const double&, const double&,
                const std::array<double, 2>&) const = 0;

    protected:
        double CrossSectionLab(bool, const double&) const noexcept;
};

class InteractionFactory {
    public:
        using TCreateMethod = std::unique_ptr<Interactions>(*)();
        TCreateMethod CreateFunc;

        static InteractionFactory& Instance();

        bool Register(const std::string&, TCreateMethod);
        std::shared_ptr<Interactions> Create(const std::string&);

    private:
        InteractionFactory() : methods() {};
        std::map<std::string, TCreateMethod> methods;
};

#define REGISTER_INTERACTION(interaction) \
    bool interaction::registered = InteractionFactory::Instance().Register(interaction::GetName(), \
            interaction::Create);

class GeantInteractions : public Interactions {
    public:
        // Initialize GeantInteractions class
        GeantInteractions(const std::string&);
        static std::unique_ptr<Interactions> Create() {
            return std::unique_ptr<GeantInteractions>(
                    new GeantInteractions("src/nuchic/data/GeantData.hdf5"));
        }

        // Destructor
        ~GeantInteractions() {};

        // Functions
        static std::string GetName() {return "GeantInteractions";}
        bool IsRegistered() const noexcept {return registered;}
        double CrossSection(const Particle&, const Particle&) const;
        ThreeVector MakeMomentum(bool, const double&, const double&,
                const std::array<double, 2>&) const;

    protected:
        static bool registered;

    private:
        // Functions
        double CrossSectionAngle(bool, const double&, const double&) const;
        void LoadData(bool, const H5::Group&);

        // Variables
        std::vector<double> m_theta, m_cdf;
        std::vector<double> m_pcmPP, m_xsecPP;
        std::vector<double> m_pcmNP, m_xsecNP;
        Interp1D m_crossSectionPP, m_crossSectionNP;
        Interp2D m_thetaDistPP, m_thetaDistNP;
};

#endif // end of include guard: INTERACTIONS_HH
