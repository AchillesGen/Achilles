#include <array>

/**
 * Get Legendre polynomials
 *
 * @param[in] x argument of Legendre Polynomials
 * @param[in] lmax maximum l
 * @param[out] P array of values of Legendre Polynomials P[l] = P_l(x)
 * @param[out] dP array of derivatives of Legendre Polynomials dP[l] = P'_l(x)
 */
//void LegendrePolynomials(double x, int lmax, std::array<double, 12> &P, std::array<double, 12> &dP){
//	lmax = std::min(lmax,11);
//
//	double x2 = x*x;
//	double x3 = x2*x;
//
//	//Starting values
//	P[0] = 1.; 
//	P[1] = x; 
//	P[2] = 1.5*x2 - 0.5;
//	P[3] = 2.5*x3 - 1.5*x;
//
//	dP[0] = 0;
//	dP[1] = 1.;
//	dP[2] = 3*x;
//	dP[3] = 7.5*x2 - 1.5;
//
//	
//	for (size_t i=4; i < lmax + 1; i++)
//	{
//		double di = static_cast <double> (i);
//		P[i] = ( (2.*di - 1.)*x*P[i-1] - (di-1.)*P[i-2]) /di;
//		dP[i] = di*P[i-1] + x*dP[i-1];
//	}
//
//}
//The former used for validation only

/**
 * Get expansion coefficients of Legendre polynomials in monomial form
 * Expansion coefficient a_l^k such that P_l(x) = \sum_k a_^l_k x^k

 * @param[in] lmax maximum l to compute (capped at 10)
 * @param[out] A where A[l][k] = a^l_k 
 */
void LegendreCoefficients(size_t lmax, std::array<std::array<double, 12>, 12> &A){
	// 	// Assuming for the moment that all elements of A are set to zero
	// Implemented up to 10
	lmax = std::min(lmax,static_cast<size_t> (10));

	A[0][0] = 1.;

	A[1][1] = 1.;

	A[2][0] = -0.5, A[2][2] = 1.5;

	A[3][1] = -1.5, A[3][3] = 2.5;

	A[4][0] = 0.125*3, A[4][2] = -3.75, A[4][4] = 4.375;

	A[5][1] = 15*0.125, A[5][3] = -8.75, A[5][5] = 7.875;

	A[6][0] = -5./16., A[6][2] = 105./16., A[6][4] = -315./16., A[6][6] = 231./16.;
	

	//Recursion for higher orders, need only the first coefficient
	A[7][1] = -35./16;
	A[9][1] = 315./128.;

	A[8][0] = 35./128.;
	A[10][0] = -63./256.;

	for (size_t n = 7 ; n < lmax + 1 ; n++)
	{
		size_t k0 = n%2 + 2;
		for (size_t k = k0 ; k <= n ; k+=2)
		{
			double dk = static_cast<double> (k);
			double dn = static_cast<double> (n);
			A[n][k] = -1.*A[n][k-2]*(dn - dk + 2)*(dn + dk - 1)/dk/(dk-1);  
		}

	}
	
}


void H_G_Coefficients(size_t n, size_t m, const std::array<std::array<double, 12>, 12> &A,  std::array<double, 26> &h, std::array<double, 26> &g)
{
	//A[n][k] are coefficients a_k such that P_n(x) = \sum_k a_k x^k
	//H = P_n (x) * P_m (x)
	//G = (1 - x^2) P_n^\prime (x) * P_m^\prime (x)
	//The output arrays h/g are coefficients such that H = \sum_k h[k]x^k , G = sum_k g[k]*x^k
	
	size_t k_n0 = n%2;
	size_t k_m0 = m%2;

	//SET arrays to zero
	//Check size of arrays
	
	//Lowest order could be done first to avoid the if-statement later

	for (size_t kn = k_n0 ; kn <= n ; kn+=2)
	{
		for (size_t km = k_m0 ; km <= m ; km+=2)
		{
			h[km+kn] += A[n][kn]*A[m][km];
			
			if (km > 0 && kn > 0)
			{
				double kmkn = static_cast<double> (km*kn);
				g[km+kn - 2] += kmkn*A[n][kn]*A[m][km];
				g[km+kn] -= -kmkn*A[n][kn]*A[m][km]; //redundant in principle, could do this in the end

			}

		}

	}
}



void H_G_Coefficients(size_t n, size_t m, const std::array<std::array<double, 12>, 12> &A,  double *h, double *g)
{
	//A[n][k] are coefficients a_k such that P_n(x) = \sum_k a_k x^k
	//H = P_n (x) * P_m (x)
	//G = (1 - x^2) P_n^\prime (x) * P_m^\prime (x)
	//The output arrays h/g are coefficients such that H = \sum_k h[k]x^k , G = sum_k g[k]*x^k
	
	size_t k_n0 = n%2;
	size_t k_m0 = m%2;

	//SET arrays to zero
	//Check size of arrays
	
	//Lowest order could be done first to avoid the if-statement later

	for (size_t kn = k_n0 ; kn <= n ; kn+=2)
	{
		for (size_t km = k_m0 ; km <= m ; km+=2)
		{
			h[km+kn] += A[n][kn]*A[m][km];
			
			if (km > 0 && kn > 0)
			{
				double kmkn = static_cast<double> (km*kn);
				g[km+kn - 2] += kmkn*A[n][kn]*A[m][km];
				g[km+kn] -= kmkn*A[n][kn]*A[m][km]; //redundant in principle, could do this in the end

			}

		}

	}
}
