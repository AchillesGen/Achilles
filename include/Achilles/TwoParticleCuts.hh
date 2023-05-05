#ifndef TWO_PARTICLE_CUT_HH
#define TWO_PARTICLE_CUT_HH

#include "Achilles/Cuts.hh"

namespace achilles {

class TwoParticleCut : public CutBase<double> {
  public:
    // Calculate cut for two input FourVectors
    virtual bool MakeCut(const FourVector &, const FourVector &) const = 0;
    TwoParticleCut(const YAML::Node &node) : CutBase(node) {}
    virtual ~TwoParticleCut() = default;
    static std::string Name() { return "Two Particle"; }
};

#define TWO_PARTICLE_CUT(CutName)                                                              \
    class CutName##Cut : public TwoParticleCut, RegistrableCut<TwoParticleCut, CutName##Cut> { \
      public:                                                                                  \
        CutName##Cut(const YAML::Node &node) : TwoParticleCut(node) {}                         \
        static std::string Name() {                                                            \
            return #CutName;                                                                   \
        }                                                                                      \
        static std::unique_ptr<TwoParticleCut> Construct(const YAML::Node &node) {             \
            return std::make_unique<CutName##Cut>(node);                                       \
        }                                                                                      \
        bool MakeCut(const FourVector &, const FourVector &) const override;                   \
    }

TWO_PARTICLE_CUT(DeltaTheta);
TWO_PARTICLE_CUT(InvariantMass);

// class DeltaThetaCut : public TwoParticleCutFactory<DeltaThetaCut> {
//     public:
//         DeltaThetaCut() = default;
//
//         static std::string Name() { return "DeltaTheta"; }
//         static std::unique_ptr<CutBase> InitializeCut(YAML::Node&);
//
//         bool MakeCut(const FourVector&, const FourVector&) override;
//
//     private:
//         cut_range m_range;
// };
//
// class InvariantMassCut : public TwoParticleCutFactory<InvariantMassCut> {
//     public:
//         InvariantMassCut() = default;
//
//         static std::string Name() { return "InvariantMass"; }
//         static std::unique_ptr<CutBase> InitializeCut(YAML::Node&);
//
//         bool MakeCut(const FourVector&, const FourVector&) override;
//
//     private:
//         cut_range m_range;
// };

} // namespace achilles

#endif
