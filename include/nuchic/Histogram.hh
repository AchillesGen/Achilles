#ifndef HISTOGRAM_HH
#define HISTOGRAM_HH

#include <string>
#include <map>
#include <mutex>
#include <stdexcept>
#include <vector>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "yaml-cpp/yaml.h"
#pragma GCC diagnostic pop

enum class HistOut { Native, YODA, ROOT };

#ifdef HAVE_YODA
#include "YODA/Histo1D.h"
#include "YODA/WriterYODA.h"
#endif

#ifdef HAVE_ROOT
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wuseless-cast"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#pragma GCC diagnostic ignored "-Wdouble-promotion"
#include "TH1.h"
#include "TFile.h"
#pragma GCC diagnostic pop
#endif


namespace nuchic {

class FourVector;

class Histogram {
    public:
        Histogram() = default;
        Histogram(const size_t&, const double&, const double&,
                  std::string, std::string path="");
        Histogram(std::vector<double>, std::string, std::string path="");
        Histogram(const Histogram&) = default;
        Histogram(Histogram&&) = default;
        Histogram& operator=(const Histogram&) = default;
        Histogram& operator=(Histogram&&) = default;
        virtual ~Histogram() = default;

        virtual void Fill(const double& x, const double& wgt=1.0);
        virtual void Scale(const double&);
        virtual void Normalize(const double& norm=1.0);
        virtual double Integral() const;
        virtual double Integral(const size_t&, const size_t&) const;
        virtual void Save() const;
        virtual void Save(const std::string&) const;

        virtual void SetName(const std::string& name_) { name = name_; }
        virtual void SetPath(const std::string& path_) { path = path_; }

        virtual std::string GetName() const { return name; }

    protected:
        std::string name, path; 
        std::vector<double> binedges;
        std::vector<double> binvals;
        std::vector<double> errors;
        size_t FindBin(const double&) const;
};

#ifdef HAVE_YODA
class YODAHistogram : public Histogram {
    public:
        YODAHistogram(const size_t&, const double&, const double&, const std::string&, const std::string& path="");
        YODAHistogram(const std::vector<double>&, const std::string&, const std::string& path="");
        void Fill(const double& x, const double& wgt=1.0);
        void Scale(const double&);
        void Normalize(const double& norm=1.0);
        double Integral() const;
        double Integral(const size_t&, const size_t&) const;
        void Save() const;
        void Save(const std::string&) const;
    private:
        YODA::Histo1D histogram;
};
#endif

#ifdef HAVE_ROOT

class ROOTHistogram : public Histogram {
    public:
        ROOTHistogram(const size_t&, const double&, const double&, const std::string&, const std::string& path="");
        ROOTHistogram(const std::vector<double>&, const std::string&, const std::string& path="");
        void Fill(const double&, const double& wgt=1.0);
        void Scale(const double&);
        void Normalize(const double& norm=1.0);
        double Integral() const;
        double Integral(const size_t&, const size_t&) const;
        void Save() const;
        void Save(const std::string&) const;
    private:
        TH1D* histogram;
};
#endif

class HistogramCollection {
    public:
        HistogramCollection(const YAML::Node&);
#if HAVE_ROOT
        ~HistogramCollection() {
            delete f;
        }
#endif

        void InitializeHists();

        bool AddHistogram(const size_t&, const double&, const double&,
                const std::string&, const std::string& path="");
        bool AddHistogram(const std::vector<double>&, const std::string&, const std::string& path="");

        bool FillHists(const std::vector<int>&, const std::vector<FourVector>&, const double&);
        bool Finalize();

        bool SaveHists() const;

    private:
#ifdef HAVE_ROOT
        TFile *f;
#endif
        HistOut outputMode;
        std::map<std::string, std::unique_ptr<Histogram>> hists;
};

}

#endif
