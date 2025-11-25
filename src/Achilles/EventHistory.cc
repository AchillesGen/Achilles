#include "Achilles/EventHistory.hh"

#include <queue>

using achilles::EventHistory;
using achilles::EventHistoryNode;

EventHistory::EventHistory(const EventHistory &other) {
    cur_idx = other.cur_idx;
    m_history.reserve(other.m_history.size());
    for(const auto &elm : other.m_history) {
        m_history.push_back(std::make_unique<EventHistoryNode>(*elm));
    }
}

EventHistory &EventHistory::operator=(const EventHistory &other) {
    if(this == &other) return *this;
    cur_idx = other.cur_idx;
    m_history.clear();
    m_history.reserve(other.m_history.size());
    for(const auto &elm : other.m_history) {
        m_history.push_back(std::make_unique<EventHistoryNode>(*elm));
    }
    return *this;
}

void EventHistory::AddVertex(ThreeVector position, const Particles &in, const Particles &out,
                             StatusCode status) {
    m_history.push_back(std::make_unique<EventHistoryNode>(cur_idx++, position, status));
    for(auto &part : in) { AddParticleIn(m_history.size() - 1, part); }
    for(auto &part : out) { AddParticleOut(m_history.size() - 1, part); }
}

void EventHistory::AddParticleIn(size_t idx, const Particle &part) {
    m_history[idx]->AddIncoming(part);
    UpdatePrevNode(part);
}

void EventHistory::AddParticleOut(size_t idx, const Particle &part) {
    m_history[idx]->AddOutgoing(part);
    UpdatePrevNode(part);
}

void EventHistory::InsertShowerVert(ThreeVector position, const Particle &org, const Particle &in,
                                    const Particle &out_org, const Particles &other) {
    // Create shower node
    m_history.push_back(
        std::make_unique<EventHistoryNode>(cur_idx++, position, StatusCode::shower));
    m_history.back()->AddIncoming(in);
    m_history.back()->AddOutgoing(out_org);
    for(auto &part : other) m_history.back()->AddOutgoing(part);

    // Insert in where original particle is located
    auto node_in = FindNodeIn(org);
    auto node_out = FindNodeOut(org);
    auto it = std::find(node_in->ParticlesIn().begin(), node_in->ParticlesIn().end(), org);
    *it = out_org;
    it = std::find(node_out->ParticlesOut().begin(), node_out->ParticlesOut().end(), org);
    *it = in;
}

std::vector<EventHistoryNode *> EventHistory::Children(size_t idx) const {
    std::vector<EventHistoryNode *> children;
    for(const auto &part : m_history[idx]->ParticlesOut()) {
        auto node = FindNodeIn(part);
        if(node) children.push_back(node);
    }
    return children;
}

std::vector<EventHistoryNode *> EventHistory::Parents(size_t idx) const {
    std::vector<EventHistoryNode *> parents;
    for(const auto &part : m_history[idx]->ParticlesIn()) {
        auto node = FindNodeOut(part);
        if(node) parents.push_back(node);
    }
    return parents;
}

std::vector<EventHistoryNode *> EventHistory::Children(EventHistoryNode *node) const {
    auto it =
        std::find_if(m_history.cbegin(), m_history.cend(),
                     [&](const std::unique_ptr<EventHistoryNode> &p) { return p.get() == node; });
    if(*it == nullptr) return {};
    return Children((*it)->Index());
}

std::vector<EventHistoryNode *> EventHistory::Parents(EventHistoryNode *node) const {
    auto it =
        std::find_if(m_history.cbegin(), m_history.cend(),
                     [&](const std::unique_ptr<EventHistoryNode> &p) { return p.get() == node; });
    if(it == m_history.end()) return {};
    return Parents((*it)->Index());
}

EventHistoryNode *EventHistory::Node(size_t idx) const {
    for(const auto &node : m_history) {
        if(node->Index() == idx) { return node.get(); }
    }
    return nullptr;
}

EventHistoryNode *EventHistory::GetUniqueNode(StatusCode status) const {
    size_t n_nodes{};
    EventHistoryNode *result = nullptr;
    for(const auto &node : m_history) {
        if(node->Status() == status) {
            result = node.get();
            n_nodes++;
        }
        if(n_nodes > 1) {
            throw std::runtime_error(
                fmt::format("EventHistory: Only one {} node is allowed!", ToString(status)));
        }
    }
    return result;
}

EventHistoryNode *EventHistory::FindNode(bool incoming, const Particle &part) const {
    for(const auto &node : m_history) {
        std::vector<Particle> particles;

        if(incoming) {
            if(node->HasIncoming(part)) return node.get();
        } else {
            if(node->HasOutgoing(part)) return node.get();
        }
    }
    return nullptr;
}

void EventHistory::UpdateStatuses(const Particles &particles) {
    for(auto &part : particles) {
        compare_momentum comp(part);
        auto node = FindNodeOut(part);
        if(node) {
            for(auto &outgoing : node->ParticlesOut()) {
                if(comp(outgoing)) outgoing.Status() = part.Status();
            }
        }
        node = FindNodeIn(part);
        if(node) {
            for(auto &incoming : node->ParticlesIn()) {
                if(comp(incoming) && incoming.Status() != ParticleStatus::target)
                    incoming.Status() = part.Status();
            }
        }
    }
}

void EventHistory::UpdatePrevNode(const Particle &part) {
    compare_momentum comp(part);
    auto node = FindNodeOut(part);
    if(!node) return;
    for(auto &outgoing : node->ParticlesOut()) {
        if(comp(outgoing)) { outgoing = part; }
    }
}

void EventHistory::WalkHistory(achilles::HistoryVisitor &visitor) const {
    // Start at primary interaction
    auto *primary = Primary();

    std::vector<size_t> visited;
    size_t current = primary->Index();
    std::queue<size_t> to_visit;

    to_visit.push(current);

    while(!to_visit.empty()) {
        current = to_visit.front();
        to_visit.pop();

        // Ensure it hasn't been visited
        if(std::find(visited.begin(), visited.end(), current) != visited.end()) continue;

        auto *node = Node(current);
        visitor.visit(node);

        for(const auto &child : Children(node)) to_visit.push(child->Index());

        for(const auto &parent : Parents(node)) { to_visit.push(parent->Index()); }

        visited.push_back(current);
    }
}
