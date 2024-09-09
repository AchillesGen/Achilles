#ifndef ACHILLES_REFHANDLER
#define ACHILLES_REFHANDLER

#include <map>
#include <string>

namespace achilles {

enum class ReferenceType {
    inspire,
    arxiv,
    doi,
};

struct Reference {
    ReferenceType type;
    std::string id;
    std::string description;
};

class ReferenceHandler {
  public:
    static ReferenceHandler &Handle() {
        static ReferenceHandler handle;
        return handle;
    }
    ~ReferenceHandler() { WriteReferences(); }

    void AddReference(Reference ref) { m_references[ref.id] = std::move(ref); }
    void WriteReferences() const;

  private:
    ReferenceHandler() = default;
    std::map<std::string, Reference> m_references;
};

} // namespace achilles

#endif
