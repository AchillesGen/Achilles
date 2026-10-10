// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef EVENT_WRITER_HH
#define EVENT_WRITER_HH

#include "Achilles/Process.hh"
#include "Achilles/Variations.hh"
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

    virtual void WriteHeader(const std::string &, const std::vector<ProcessGroup> &) = 0;
    virtual void Write(const Event &) = 0;
    /// Variation groups and their weights, must be set before writing the header
    void SetVariations(std::vector<VariationGroupInfo> info) {
        m_variation_info = std::move(info);
        m_variation_names.clear();
        for(const auto &group : m_variation_info)
            m_variation_names.insert(m_variation_names.end(), group.weight_names.begin(),
                                     group.weight_names.end());
    }

  protected:
    static std::ostream *InitializeStream(const std::string &, bool);
    std::vector<VariationGroupInfo> m_variation_info{};
    std::vector<std::string> m_variation_names{};
};

class AchillesWriter : public EventWriter {
  public:
    AchillesWriter(const std::string &filename, bool zip = true)
        : m_out{EventWriter::InitializeStream(filename, zip)}, toFile{true}, zipped{zip} {}
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

    void WriteHeader(const std::string &, const std::vector<ProcessGroup> &) override;
    void Write(const Event &) override;

  private:
    std::ostream *m_out;
    bool toFile{false};
    bool zipped{true};
    size_t nEvents{0};
};

} // namespace achilles

#endif
