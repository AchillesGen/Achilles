#ifndef F_POLYNOMIAL_H
#define F_POLYNOMIAL_H

class f_Polynomial {
  public:
    f_Polynomial();
    f_Polynomial(const double *coeffs, int _order) { SetCoeffs(coeffs, _order); }
    f_Polynomial(const f_Polynomial &f) { SetCoeffs(f.Coeffs, f.order); }

    void SetCoeffs(const double *coeffs, int _order);
    double eval_f(double x);
    double eval_df(double x);
    double eval_If(double x);

    //	    void print();

  private:
    const static int max_k = 50;
    int order = 0;
    double Coeffs[max_k] = {
        0.}; // Could do complex if neccesary, but all coeffs should be real in principle
};

#endif
