#ifndef NUCLEON_PDF_HH
#define NUCLEON_PDF_HH

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wold-style-cast"
#include "LHAPDF/LHAPDF.h"
#pragma GCC diagnostic pop
#include <string>

namespace nuchic {

class PID;

class NucleonPDF {
    public:
        NucleonPDF(const std::string&, int=0);
        ~NucleonPDF() {
            if(m_pdf) {
                delete m_pdf;
                m_pdf = nullptr;
            }
            for(auto &pdf : m_pdfs) {
                delete pdf;
            }
        }

        // Get set information
        size_t NSets() const { return m_set.size(); }

        // Access raw pdfs
        double xfxQ2(const nuchic::PID&, double, double) const;
        double xfxQ2(int, double, double) const;
        double xfxQ2(size_t, const nuchic::PID&, double, double) const;
        double xfxQ2(size_t, int, double, double) const;

        // Obtain uncertainty
        LHAPDF::PDFUncertainty Uncertainty(const std::vector<double>&);

        // Calculate structure functions
        double F2(const nuchic::PID&, double, double, size_t=0) const;
        double F1(const nuchic::PID&, double, double, size_t=0) const;
        std::vector<double> F2All(const nuchic::PID&, double, double, size_t=0) const;
        std::vector<double> F1All(const nuchic::PID&, double, double, size_t=0) const;

    private:
        double C2q1xFPlus(int, double, double, double) const;
        double C2q1xF(int, double, double) const;
        double C2g1xF(double, double, double) const;
        double C1q1xF(int, double, double, double) const;
        double C1g1xF(double, double, double) const;
        double P1qqPlus(int, double, double, double) const;
        double P1qq(int, double, double) const;
        double P1qg(double, double, double) const;

        std::string m_name;
        int m_iset;
        LHAPDF::PDFSet m_set;
        const LHAPDF::PDF* m_pdf{nullptr};
        const std::vector<LHAPDF::PDF*> m_pdfs;
        static constexpr double CF = 4.0/3.0;
        static constexpr double TR = 1.0/2.0;
        static constexpr double NC = 3;
        static constexpr double NF = 4;
};

}

#endif
