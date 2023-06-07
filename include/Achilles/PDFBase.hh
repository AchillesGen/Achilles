#ifndef PDFBASE_HH
#define PDFBASE_HH

namespace achilles {

class PDFBase {
    public:
        virtual ~PDFBase() = default;
        virtual double operator()(double, double) const = 0;
};

}


#endif
