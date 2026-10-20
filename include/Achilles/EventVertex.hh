#ifndef EVENT_VERTEX_HH
#define EVENT_VERTEX_HH

#include "Achilles/ThreeVector.hh"

namespace achilles {

class EventVertex {
  public:
    enum class StatusCode : int {
        propagation = -2,
        cascade = -1,
        primary = 1,
        beam = 2,
        target = 3,
        decay = 4,
        shower = 5,
    };

    EventVertex(size_t idx,StatusCode status=StatusCode::cascade,ThreeVector position={})
        : m_idx{idx}, m_status{status}, m_position{position} {}
    EventVertex(const EventVertex &) = default;

    // Status functions
    StatusCode &Status() { return m_status; }
    const StatusCode &Status() const { return m_status; }
    bool IsCascade() const { return m_status == StatusCode::cascade; }
    bool IsPropagation() const { return m_status == StatusCode::propagation; }
    bool IsDecay() const { return m_status == StatusCode::decay; }
    bool IsInternal() const { return IsCascade() || IsPropagation(); }
    bool IsExternal() const { return !IsInternal(); }

    // Graph accessors
    size_t Index() const { return m_idx; }
    ThreeVector& Position() { return m_position; }
  private:
    size_t m_idx{};
    ThreeVector m_position{};
    StatusCode m_status;
};

inline std::string ToString(EventVertex::StatusCode code) {
    using StatusCode = EventVertex::StatusCode;
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

} // namespace achilles

#endif