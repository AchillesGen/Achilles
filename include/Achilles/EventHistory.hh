#ifndef EVENT_HISTORY_HH
#define EVENT_HISTORY_HH

#include <memory>

#include "Achilles/Particle.hh"
#include "Achilles/ParticleInfo.hh"

namespace achilles {

class Nucleus;
class NuclearRemnant;

struct compare_momentum : public std::unary_function<Particle, bool> {
    explicit compare_momentum(const Particle &particle, double _eps = 1e-10) : self(particle), eps(_eps) {}
    bool operator()(const Particle &other) {
        return self.Momentum().Approx(other.Momentum(), eps) && self.ID() == other.ID();
    }
    Particle self;
    double eps;
};

class EventHistoryNode {
    public:
        enum class StatusCode {
            propagation = -2,
            cascade = -1,
            primary = 0,
            beam = 1,
            target = 2,
            decay = 3,
            shower = 4,
        };

        EventHistoryNode(size_t idx, StatusCode status = StatusCode::cascade) : m_idx{idx}, m_status{status} {}
        EventHistoryNode(size_t idx, ThreeVector position, StatusCode status = StatusCode::cascade) 
            : m_idx{idx}, m_position{position}, m_status{status} {}
        void AddIncoming(const Particle &part) { m_particles_in.push_back(part); }
        void AddOutgoing(const Particle &part) { m_particles_out.push_back(part); }

        // Status functions
        StatusCode& Status() { return m_status; }
        const StatusCode& Status() const { return m_status; }
        bool IsCascade() const { return m_status == StatusCode::cascade; }
        bool IsPropagation() const { return m_status == StatusCode::propagation; }
        bool IsDecay() const { return m_status == StatusCode::decay; }
        bool IsInternal() const { return IsCascade() || IsPropagation(); }
        bool IsExternal() const { return !IsInternal(); }

        // Graph accessors
        size_t Index() const { return m_idx; }
        bool HasIncoming(const Particle &part) { 
            return std::find_if(m_particles_in.begin(), m_particles_in.end(),
                                compare_momentum(part)) != m_particles_in.end();
        }
        bool HasOutgoing(const Particle &part) {
            return std::find_if(m_particles_out.begin(), m_particles_out.end(),
                                compare_momentum(part)) != m_particles_out.end();
        }
        const std::vector<Particle>& ParticlesIn() const { return m_particles_in; }
        std::vector<Particle>& ParticlesIn() { return m_particles_in; }
        const std::vector<Particle>& ParticlesOut() const { return m_particles_out; }
        std::vector<Particle>& ParticlesOut() { return m_particles_out; }
        const ThreeVector& Position() const { return m_position; }

    private:
        std::vector<Particle> m_particles_in{}, m_particles_out{};
        size_t m_idx;
        ThreeVector m_position{};
        StatusCode m_status;
};

inline std::string ToString(EventHistoryNode::StatusCode code) {
    using StatusCode = EventHistoryNode::StatusCode;
    switch(code) {
        case StatusCode::propagation:
            return "propagation";
        case StatusCode::cascade:
            return "cascade";
        case StatusCode::primary:
            return "primary";
        case StatusCode::beam:
            return "beam";
        case StatusCode::target:
            return "target";
        case StatusCode::decay:
            return "decay";
        case StatusCode::shower:
            return "shower";
    }
    return "UNKOWN";
}

class HistoryVisitor {
    public:
        virtual ~HistoryVisitor() = default;
        virtual void visit(EventHistoryNode*) = 0;
};

class EventHistory {
    public:
        using StatusCode = EventHistoryNode::StatusCode;
        EventHistory() = default;
        void AddVertex(ThreeVector position, const std::vector<Particle> &in = {},
                       const std::vector<Particle> &out = {}, StatusCode status = StatusCode::cascade);
        void AddParticleIn(size_t idx, const Particle &part);
        void AddParticleOut(size_t idx, const Particle &part);
        void InsertShowerVert(ThreeVector position, const Particle &org, const Particle &in,
                              const Particle &out_org, const std::vector<Particle> &other);

        // Accessors
        EventHistoryNode* Node(size_t idx) const;
        std::vector<EventHistoryNode*> Children(size_t idx) const;
        std::vector<EventHistoryNode*> Parents(size_t idx) const;
        std::vector<EventHistoryNode*> Children(EventHistoryNode*) const;
        std::vector<EventHistoryNode*> Parents(EventHistoryNode*) const;
        EventHistoryNode* FindNodeIn(const Particle &part) const { return FindNode(true, part); }
        EventHistoryNode* FindNodeOut(const Particle &part) const { return FindNode(false, part); }
        EventHistoryNode* Primary() const { return GetUniqueNode(StatusCode::primary); }
        EventHistoryNode* Beam() const { return GetUniqueNode(StatusCode::beam); }
        EventHistoryNode* Target() const { return GetUniqueNode(StatusCode::target); }

        // Transversal
        void WalkHistory(HistoryVisitor&) const;

        // Information
        size_t size() { return m_history.size(); }

    private:
        EventHistoryNode* FindNode(bool, const Particle&) const;
        EventHistoryNode* GetUniqueNode(StatusCode) const;
        std::vector<std::unique_ptr<EventHistoryNode>> m_history{};
        size_t cur_idx{};
};

struct PrintVisitor : HistoryVisitor {
    std::string data;
    void visit(achilles::EventHistoryNode *node) {
        data += fmt::format("Node({}, {{{}}} -> {{{}}})\n",
                            achilles::ToString(node -> Status()),
                            fmt::join(node -> ParticlesIn().begin(),
                                      node -> ParticlesIn().end(), ", "),
                            fmt::join(node -> ParticlesOut().begin(),
                                      node -> ParticlesOut().end(), ", "));
    }
};

struct PIDLocator : HistoryVisitor {
    std::vector<Particle> particles;
    PID pid;
    int direction;
    PIDLocator(PID _pid, int dir) : pid(_pid), direction(dir) {}
    void visit(achilles::EventHistoryNode *node) {
        if(direction >= 0) {
            for(auto &particle : node -> ParticlesOut()) {
                if(particle.ID() == pid) particles.push_back(particle);
            }
        } else {
            for(auto &particle : node -> ParticlesIn()) {
                if(particle.ID() == pid) particles.push_back(particle);
            }
        }
    }
};

}

#endif
