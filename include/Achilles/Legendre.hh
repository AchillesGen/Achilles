#ifndef LEGENDRE_HH
#define LEGENDRE_HH

#include <array>

//void LegendrePolynomials(double x, int lmax, std::array<double, 12> &P, std::array<double, 12> &dP);

void LegendreCoefficients(size_t lmax, std::array<std::array<double, 12>, 12> &A);


void H_G_Coefficients(size_t n, size_t m, const std::array<std::array<double, 12>, 12> &A,  std::array<double, 26> &h, std::array<double, 26> &g);


void H_G_Coefficients(size_t n, size_t m, const std::array<std::array<double, 12>, 12> &A,  double *h, double *g);

#endif