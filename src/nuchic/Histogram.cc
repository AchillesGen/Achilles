#include "nuchic/Histogram.hh"
#include "fmt/format.h"

#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <iostream>

namespace nuchic {

Histogram::Histogram(const size_t& nbins, const double& lower, const double& upper,
                     std::string name_, std::string path_) :
                    name{std::move(name_)}, path{std::move(path_)} {
    binvals = std::vector<double>(nbins, 0);
    errors = std::vector<double>(nbins, 0);
    double binsize = (upper-lower)/static_cast<double>(nbins);
    binedges.resize(nbins+1);
    for(size_t i = 0; i <= nbins; ++i) {
        binedges[i] = lower+static_cast<double>(i)*binsize;
    }
}

Histogram::Histogram(std::vector<double> binedges_, std::string name_,
                     std::string path_) : 
                    name{std::move(name_)}, path{std::move(path_)}, binedges{std::move(binedges_)} {
    binvals = std::vector<double>(binedges.size()-1, 0);
    errors = std::vector<double>(binedges.size()-1, 0);
    for(size_t i = 0; i < binvals.size(); ++i) {
        if(binedges[i]>binedges[i+1]) throw std::runtime_error("Bin edges must be in increasing order!");
        if(binedges[i]==binedges[i+1]) throw std::runtime_error("Bin edges must not be identical!");
    }
}

size_t Histogram::FindBin(const double& x) const {
    auto it = std::lower_bound(binedges.begin(),binedges.end(), x);
    if(it == binedges.end()) return static_cast<size_t>(-1);

    return static_cast<size_t>(std::distance(binedges.begin(),it));
}

void Histogram::Fill(const double& x, const double& wgt) {
    const size_t loc = FindBin(x);
    nentries++;
    if(loc != static_cast<size_t>(-1) && loc != 0) {
        binvals[loc-1] += wgt/(binedges[loc]-binedges[loc-1]);
        errors[loc-1] += pow(wgt/(binedges[loc]-binedges[loc-1]), 2);
    }
}

void Histogram::Scale(const double& scale) {
    for(double & binval : binvals) {
        binval *= scale;
    }
}

void Histogram::Normalize(const double& norm) {
    Scale(norm/Integral());
}

double Histogram::Integral() const {
    return Integral(0, binedges.size());
}

double Histogram::Integral(const size_t& lower, const size_t& upper) const {
    double result = 0;
    double sign = 1;
    if(lower > upper)
        sign = -1;

    if(lower > binedges.size())
        throw std::runtime_error("Invalid range for histogram integration");
    if(upper > binedges.size())
        throw std::runtime_error("Invalid range for histogram integration");

    for(size_t i = lower; i < upper; ++i) {
        const auto mean = binvals[i];///static_cast<double>(nentries);
        result += mean*(binedges[i+1]-binedges[i]);
    }

    return sign*result;
}

void Histogram::Save(std::ostream *out) const {
    *out << name << std::endl;
    *out << fmt::format("{:^15} {:^15} {:^15} {:^15}\n",
                       "lower edge", "upper edge", "value", "error");
    for(size_t i = 0; i < binvals.size(); ++i) {
        *out << fmt::format("{:< 15.6e} {:< 15.6e} {:< 15.6e} {:< 15.6e}\n",
                           binedges[i], binedges[i+1], binvals[i], std::sqrt(errors[i]));
    }
}

void Histogram::Save(const std::string& filename) const {
    std::ofstream out(path+filename+".txt");
    Save(&out);
}

#ifdef HAVE_YODA

YODAHistogram::YODAHistogram(const size_t& nbins, const double& lower, const double& upper, const std::string& name_, const std::string& path_) {
    name = name_;
    path = path_;
    histogram = YODA::Histo1D(nbins,lower,upper,path,name);
}

YODAHistogram::YODAHistogram(const std::vector<double>& binedges, const std::string& name_, const std::string& path_) {
    name = name_;
    path = path_;
    histogram = YODA::Histo1D(binedges,name,path);
}

void YODAHistogram::Fill(const double& x, const double& wgt) {
    histogram.fill(x,wgt);
}

void YODAHistogram::Scale(const double& scale) {
    histogram.scaleW(scale);
}

void YODAHistogram::Normalize(const double& norm) {
    histogram.normalize(norm);
}

double YODAHistogram::Integral() const {
    return histogram.integral();
}

double YODAHistogram::Integral(const size_t& low, const size_t& high) const {
    return histogram.integralRange(low,high);
}

void YODAHistogram::Save() const {
    std::string filename = path + name +".yoda";
    YODA::WriterYODA::write(filename,histogram);
}

void YODAHistogram::Save(const std::string& filename) const {
    std::string outfile = path + filename +".yoda";
    YODA::WriterYODA::write(outfile,histogram);
}

#endif

#ifdef HAVE_ROOT

ROOTHistogram::ROOTHistogram(const size_t& nbins, const double& lower, const double& upper, const std::string& _name, const std::string&) {
    histogram = new TH1D(_name.c_str(), "", static_cast<int>(nbins), lower, upper);
}

ROOTHistogram::ROOTHistogram(const std::vector<double>& _binedges, const std::string& _name, const std::string&) {
    const size_t nbins = _binedges.size()-1;
    std::vector<double> bins(nbins);
    for(size_t i = 0; i <= nbins; ++i) bins[i] = _binedges[i];
    histogram = new TH1D(_name.c_str(), "", static_cast<int>(nbins), bins.data());
}

void ROOTHistogram::Fill(const double& x, const double& wgt) {
    double bwidth = histogram -> GetBinWidth(histogram -> FindFixBin(x));
    histogram -> Fill(x,wgt/bwidth);
}

void ROOTHistogram::Scale(const double& scale) {
    histogram -> Scale(scale);
}

void ROOTHistogram::Normalize(const double& norm) {
    histogram -> Scale(norm/Integral());
}

double ROOTHistogram::Integral() const {
    return histogram -> Integral();
}

double ROOTHistogram::Integral(const size_t& low, const size_t& high) const {
    return histogram -> Integral(static_cast<int>(low), static_cast<int>(high));
}

void ROOTHistogram::Save() const {
}

void ROOTHistogram::Save(const std::string&) const {
}

#endif

/*
HistogramCollection::HistogramCollection(IO::Settings *settings) {
    std::string mode = settings -> GetSettingString("HistogramMode");
    if(mode == "BuiltIn") outputMode = HistOut::Native;
    else if(mode == "YODA") outputMode = HistOut::YODA;
    else if(mode == "ROOT") outputMode = HistOut::ROOT;
    else throw std::runtime_error("Output mode " + mode + " is not supported");

#ifdef HAVE_ROOT
    if(outputMode == HistOut::ROOT) {
        f = new TFile("resbos.root","RECREATE");
    }
#endif

    InitializeHists();
}

bool HistogramCollection::AddHistogram(const size_t& nbins, const double& lower,
        const double& upper, const std::string& name, const std::string& path) {
    Histogram *hist;
    switch(outputMode) {
        case HistOut::Native:
            hist = new Histogram(nbins,lower,upper,name,path); 
            break;

        case HistOut::YODA:
#ifdef HAVE_YODA
            hist = new YODAHistogram(nbins,lower,upper,name,path);
            break;
#else
            throw std::runtime_error("ResBos2 not compiled with YODA histogram support");
#endif

        case HistOut::ROOT:
#ifdef HAVE_ROOT
            hist = new ROOTHistogram(nbins,lower,upper,name,path);
            break;
#else
            throw std::runtime_error("ResBos2 not compiled with ROOT histogram support");
#endif

        default:
            throw std::runtime_error("Invalid histogramming mode");
    }
    hists[name] = hist;
    return true;
}

bool HistogramCollection::AddHistogram(const std::vector<double>& bins,
        const std::string& name, const std::string& path) {
    Histogram *hist = nullptr;
    switch(outputMode) {
        case HistOut::Native:
            hist = new Histogram(bins,name,path); 
            break;

        case HistOut::YODA:
#ifdef HAVE_YODA
            hist = new YODAHistogram(bins,name,path);
            break;
#else
            throw std::runtime_error("ResBos2 not compiled with YODA histogram support");
#endif

        case HistOut::ROOT:
#ifdef HAVE_ROOT
            hist = new ROOTHistogram(bins,name,path);
            break;
#else
            throw std::runtime_error("ResBos2 not compiled with ROOT histogram support");
#endif

        default:
            throw std::runtime_error("Invalid histogramming mode");
    }
    hists[name] = hist;
    return true;
}

bool HistogramCollection::SaveHists() const {
#ifdef HAVE_ROOT
    if(outputMode == HistOut::ROOT) {
        f -> Write();
        f -> Close();
    } else {
#endif
    for(auto hist : hists) {
        hist.second -> Save(hist.second->GetName());
    }
#ifdef HAVE_ROOT
    }
#endif
    return true;
}
*/

}
