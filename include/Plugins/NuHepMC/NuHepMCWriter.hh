// SPDX-FileCopyrightText: 2018-2026 Achilles Developers
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef NUHEPMC_WRITER
#define NUHEPMC_WRITER

#include "Achilles/EventWriter.hh"
#include "Achilles/Statistics.hh"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wconversion"
#include "HepMC3/Writer.h"
#pragma GCC diagnostic pop

namespace achilles {

class ProcessGroup;

class NuHepMCWriter : public EventWriter {
  public:
    NuHepMCWriter(const string &filename, bool) : outfilename{filename} {}
    ~NuHepMCWriter() override = default;

    void WriteHeader(const std::string &, const std::vector<ProcessGroup> &) override;
    void Write(const Event &) override;
    void SetWeightNormalization(double norm) override { m_weight_norm = norm; }

  private:
    std::shared_ptr<HepMC3::Writer> file;
    std::string outfilename;
    achilles::StatsData results;
    std::vector<achilles::StatsData> variation_results;
    double m_weight_norm{1};
    static constexpr std::array<int, 3> version{0, 1, 0};
};

} // namespace achilles

#endif
