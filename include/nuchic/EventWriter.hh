#ifndef EVENT_WRITER_HH
#define EVENT_WRITER_HH

#include <ostream>
#include <fstream>
#include <string>
#include <vector>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "gzstream/gzstream.h"
#pragma GCC diagnostic pop

namespace nuchic {

class Event;

class EventWriter {
    public:
        EventWriter() = default;
        EventWriter(const EventWriter&) = default;
        EventWriter(EventWriter&&) = default;
        EventWriter& operator=(const EventWriter&) = default;
        EventWriter& operator=(EventWriter&&) = default;
        virtual ~EventWriter() = default;

        virtual void WriteHeader(const std::string&) = 0;
        virtual void Write(const Event&) = 0;
};

class NuchicWriter : public EventWriter {
    public:
        NuchicWriter(const std::string&, bool=true);
        NuchicWriter(std::ostream *out) : m_out{out} {}
        NuchicWriter(const NuchicWriter&) = default;
        NuchicWriter(NuchicWriter&&) = default;
        NuchicWriter& operator=(const NuchicWriter&) = default;
        NuchicWriter& operator=(NuchicWriter&&) = default;
        ~NuchicWriter() override {
            if(toFile) {
                if(zipped) {
                    dynamic_cast<ogzstream*>(m_out) -> close();
                } else {
                    dynamic_cast<std::ofstream*>(m_out) -> close();
                }
                delete m_out;
            }
            m_out = nullptr;
        }

        void WriteHeader(const std::string&) override;
        void Write(const Event&) override;

    private:
        bool toFile{false}, zipped{true};
        size_t nEvents{0};
        std::ostream *m_out; 
};

}

#endif
