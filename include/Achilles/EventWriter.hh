#ifndef EVENT_WRITER_HH
#define EVENT_WRITER_HH

#include <fstream>
#include <ostream>
#include <string>
#include <vector>

#if GZIP
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "gzstream/gzstream.h"
#pragma GCC diagnostic pop
#endif

namespace achilles {

class Event;

class EventWriter {
  public:
    EventWriter() = default;
    EventWriter(const EventWriter &) = default;
    EventWriter(EventWriter &&) = default;
    EventWriter &operator=(const EventWriter &) = default;
    EventWriter &operator=(EventWriter &&) = default;
    virtual ~EventWriter() = default;

    virtual void WriteHeader(const std::string &) = 0;
    virtual void Write(const Event &) = 0;
};

class AchillesWriter : public EventWriter {
  public:
    AchillesWriter(const std::string &, bool = true);
    AchillesWriter(std::ostream *out) : m_out{out} {}
    AchillesWriter(const AchillesWriter &) = default;
    AchillesWriter(AchillesWriter &&) = default;
    AchillesWriter &operator=(const AchillesWriter &) = default;
    AchillesWriter &operator=(AchillesWriter &&) = default;
    ~AchillesWriter() override {
        if(toFile) {
#ifdef GZIP
            if(zipped) {
                dynamic_cast<ogzstream *>(m_out)->close();
            } else {
                dynamic_cast<std::ofstream *>(m_out)->close();
            }
#else
            dynamic_cast<std::ofstream *>(m_out)->close();
#endif
            delete m_out;
        }
        m_out = nullptr;
    }

    void WriteHeader(const std::string &) override;
    void Write(const Event &) override;

  private:
    bool toFile{false};
    bool zipped{true};
    size_t nEvents{0};
    std::ostream *m_out;
};

} // namespace achilles

#endif
