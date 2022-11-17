#ifndef EVENT_HISTORY_HH
#define EVENT_HISTORY_HH

#include <memory>

#include "Achilles/Particle.hh"
#include "Achilles/ParticleInfo.hh"

namespace achilles {

class Nucleus;
class NuclearRemnant;

class EventHistoryNode {
    public:
        enum class StatusCode {
            propagation = -2,
            cascade = -1,
            primary = 0,
            beam = 1,
            target = 2,
            external = 3,
            decay = 4,
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
        bool IsInternal() const { return IsCascade() || IsPropagation(); }
        bool IsExternal() const { return !IsInternal(); }

        // Graph accessors
        size_t Index() const { return m_idx; }
        bool HasIncoming(const Particle &part) { 
            return std::find(m_particles_in.begin(), m_particles_in.end(), part) != m_particles_in.end();
        }
        bool HasOutgoing(const Particle &part) {
            return std::find(m_particles_out.begin(), m_particles_out.end(), part) != m_particles_out.end();
        }
        const std::vector<Particle>& ParticlesIn() const { return m_particles_in; }
        const std::vector<Particle>& ParticlesOut() const { return m_particles_out; }

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
        case StatusCode::external:
            return "external";
        case StatusCode::decay:
            return "decay";
    }
    return "UNKOWN";
}

class EventHistory {
    public:
        using StatusCode = EventHistoryNode::StatusCode;
        EventHistory() = default;
        void AddVertex(ThreeVector position, const std::vector<Particle> &in = {},
                       const std::vector<Particle> &out = {}, StatusCode status = StatusCode::cascade) {
            m_history.push_back(std::make_unique<EventHistoryNode>(cur_idx++, position, status));
            for(const auto &part : in) {
                m_history.back() -> AddIncoming(part);
            }
            for(const auto &part : out) {
                m_history.back() -> AddOutgoing(part);
            }
        }
        EventHistoryNode* Node(size_t idx) const;
        void AddParticleIn(size_t idx, const Particle &part);
        void AddParticleOut(size_t idx, const Particle &part);
        std::vector<EventHistoryNode*> Children(size_t idx) const;
        std::vector<EventHistoryNode*> Parents(size_t idx) const;
        std::vector<EventHistoryNode*> Children(EventHistoryNode*) const;
        std::vector<EventHistoryNode*> Parents(EventHistoryNode*) const;
        EventHistoryNode* FindNodeIn(const Particle &part) const { return FindNode(true, part); }
        EventHistoryNode* FindNodeOut(const Particle &part) const { return FindNode(false, part); }
        size_t size() { return m_history.size(); }
        EventHistoryNode* Primary() const { return GetUniqueNode(StatusCode::primary); }
        EventHistoryNode* Beam() const { return GetUniqueNode(StatusCode::beam); }
        EventHistoryNode* Target() const { return GetUniqueNode(StatusCode::target); }
 
    private:
        EventHistoryNode* FindNode(bool, const Particle&) const;
        EventHistoryNode* GetUniqueNode(StatusCode) const;
        std::vector<std::unique_ptr<EventHistoryNode>> m_history{};
        size_t cur_idx{};
};

}

#endif
