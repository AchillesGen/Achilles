#ifndef EVENT_HISTORY_HH
#define EVENT_HISTORY_HH

#include "Achilles/Particle.hh"

namespace achilles {

class EventHistory {
    public:
        EventHistory() = default;

        size_t Length() const { return m_vertices.size(); }
        void AddVertex();
        void AddIncoming(size_t);
        void AddOutgoing(size_t);
        void ConnectVertex(size_t, size_t, size_t, size_t);

    private:
        using Edge = std::vector<std::pair<Particle, Particle>>;
        struct Vertex {
            std::vector<Particle> in, out;
            ThreeVector location;
        };

        std::vector<Vertex> m_vertices;
        std::vector<Edge> m_edges;
};

}

#endif
