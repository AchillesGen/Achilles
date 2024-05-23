#pragma once

#include "Achilles/Cascade.hh"
#include "Achilles/Histogram.hh"
#include "yaml-cpp/yaml.h"

#include <string>

namespace achilles {

class Nucleus;

enum class CascadeMode {
    CrossSection,
    MeanFreePath,
    Transparency,
    CrossSectionMFP,
    TransparencyMFP,
    TransparencyExternal
};

} // namespace achilles

namespace YAML {

template <> struct convert<achilles::CascadeMode> {
    static bool decode(const YAML::Node &node, achilles::CascadeMode &mode) {
        std::string name = node.as<std::string>();
        if(name == "CrossSection")
            mode = achilles::CascadeMode::CrossSection;
        else if(name == "MeanFreePath")
            mode = achilles::CascadeMode::MeanFreePath;
        else if(name == "Transparency")
            mode = achilles::CascadeMode::Transparency;
        else if(name == "CrossSectionMFP")
            mode = achilles::CascadeMode::CrossSectionMFP;
        else if(name == "TransparencyMFP")
            mode = achilles::CascadeMode::TransparencyMFP;
	else if(name == "TransparencyExternal")
            mode = achilles::CascadeMode::TransparencyExternal;
        else
            return false;

        return true;
    }
};

} // namespace YAML

namespace achilles::CascadeTest {

void RunCascade(const std::string &runcard);

class RunMode {
  public:
    RunMode(std::shared_ptr<Nucleus> nuc, Cascade cascade)
        : m_nuc{nuc}, m_cascade{std::move(cascade)} {}
    virtual ~RunMode() = default;
    virtual void GenerateEvent(double) = 0;
    virtual void PrintResults(std::ofstream &) const = 0;
    virtual void Reset() = 0;

  protected:
    std::shared_ptr<Nucleus> m_nuc;
    Cascade m_cascade;
};

class CalcCrossSection : public RunMode {
  public:
    CalcCrossSection(int pid, std::shared_ptr<Nucleus> nuc, Cascade cascade, double radius = 10)
        : RunMode(nuc, std::move(cascade)), m_radius{std::move(radius)}, m_pid{pid} {}
    void GenerateEvent(double mom) override;
    void PrintResults(std::ofstream &out) const override;
    void Reset() override {
        nevents = 0;
        nhits = 0;
    }

  private:
    double m_radius;
    int m_pid;
    double nhits{}, nevents{};
};

class CalcCrossSectionMFP : public RunMode {
  public:
    CalcCrossSectionMFP(int pid, std::shared_ptr<Nucleus> nuc, Cascade cascade, double radius = 10)
        : RunMode(nuc, std::move(cascade)), m_radius{std::move(radius)}, m_pid{pid} {}
    void GenerateEvent(double mom) override;
    void PrintResults(std::ofstream &out) const override;
    void Reset() override {
        nevents = 0;
        nhits = 0;
    }

  private:
    double m_radius;
    int m_pid;
    double nhits{}, nevents{};
};

class CalcMeanFreePath : public RunMode {
  public:
    CalcMeanFreePath(int pid, std::shared_ptr<Nucleus> nuc, Cascade cascade);
    void GenerateEvent(double kick_mom) override;
    void PrintResults(std::ofstream &out) const override;

    // TODO: Implement this!
    void Reset() override {}

  private:
    int m_pid;
    Histogram m_hist;
};

class CalcTransparency : public RunMode {
  public:
    CalcTransparency(std::shared_ptr<Nucleus> nuc, Cascade cascade)
        : RunMode(nuc, std::move(cascade)) {}

    void GenerateEvent(double kick_mom) override;
    void PrintResults(std::ofstream &out) const override;
    void Reset() override {
        nevents = 0;
        ninteract = 0;
        ncaptured = 0;
    }

  private:
    double nevents{};
    double ninteract{};
    double distance{};
    double ncaptured{};
};

class CalcTransparencyMFP : public RunMode {
  public:
    CalcTransparencyMFP(std::shared_ptr<Nucleus> nuc, Cascade cascade)
        : RunMode(nuc, std::move(cascade)) {}

    void GenerateEvent(double kick_mom) override;
    void PrintResults(std::ofstream &out) const override;
    void Reset() override {
        nevents = 0;
        ninteract = 0;
        distance = 0;
        ncaptured = 0;
    }

  private:
    double nevents{};
    double ninteract{};
    double distance{};
    double ncaptured{};
};


class CalcTransparencyExternal : public RunMode {
  public:
    CalcTransparencyExternal(int pid, std::shared_ptr<Nucleus> nuc, Cascade cascade)
        : RunMode(nuc, std::move(cascade)), m_pid{pid} {}
    void GenerateEvent(double mom) override;
    void PrintResults(std::ofstream &out) const override;
    void Reset() override {
        nevents = 0;
        ninteract = 0;
        distance = 0;
        ncaptured = 0;
    }

  private:
    int m_pid;
    double nevents{};
    double ninteract{};
    double distance{};
    double ncaptured{};
};

} // namespace achilles::CascadeTest
