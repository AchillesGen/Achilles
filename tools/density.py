import numpy as np
import argparse
from scipy.special import jv


FB_6_6A = np.array([8.0, 0.15721e-1, 0.38732e-1, 0.36808e-1, 0.14671e-1,
                    -0.43277e-2, -0.97752e-2, -0.68908e-2, -0.27631e-2,
                    -0.63568e-3, 0.71809e-5, 0.18441e-3, 0.75066e-4,
                    0.51069e-4, 0.14308e-4, 0.23170e-5, 0.68465e-6, 0])

FB_6_6B = np.array([8.0, 0.15737e-1, 0.38897e-1, 0.37085e-1, 0.14795e-1,
                    -0.44831e-2, -0.10057e-1, -0.68696e-2, -0.28813e-2,
                    -0.77229e-3, 0.66908e-4,  0.10636e-3, -0.36864e-4,
                    -0.50135e-5, 0.94550e-5, -0.47687e-5, 0, 0])

FB_12_18 = np.array([9.0, 0.30451e-1, 0.55337e-1, 0.20203e-1, -0.16765e-1,
                     -0.13578e-1, -0.43204e-4, 0.91988e-3, -0.41205e-3,
                     0.11971e-3, -0.19801e-4, -0.43204e-5, 0.61205e-5,
                     -0.37803e-5, 0.18001e-5, -0.77407e-6, 0, 0])

SOG_6_6 = np.array([2.4696, 1.20, 0.0, 0.016690, 0.4, 0.050325, 1.0, 0.128621,
                    1.3, 0.180515, 1.7, 0.219097, 2.3, 0.278416, 2.7, 0.058779,
                    3.5, 0.057817, 4.3, 0.007739, 5.4, 0.002001, 6.7, 0.000007,
                    0, 0])


def woods_saxon(r, A):
    R = (2.745e-4*A+1.063)*A**(1/3)
    D = 0.510 + 1.63e-4*A
    return 1/(1+np.exp((r-R)/D))


def fb_model(r, params):
    results = np.zeros_like(r)
    x = np.pi*r/params[0]
    for i, param in enumerate(params[1:]):
        results += param*jv(0, i*x)
    return np.maximum(0, results)


def sog_model(r, params):
    fact = 2*np.pi*np.sqrt(np.pi)
    g = params[0]/np.sqrt(1.5)
    coeff = fact*g**3
    results = np.zeros_like(r)
    for i in range(2, 25, 2):
        ri = params[i]
        qi = params[i+1]
        ai = qi/(coeff*(1+2*(ri/g)**2))

        results += ai*(np.exp(-((r-ri)/g)**2)+np.exp(-((r+ri)/g)**2))

    return np.maximum(0, results)


def main():
    parser = argparse.ArgumentParser(description="Density profile calculation")
    parser.add_argument('nucleons', metavar='A', type=int,
                        help='Number of nucleons in the nucleus')
    parser.add_argument('--r0', type=float,
                        help='Minimum radius to consider', default=0)
    parser.add_argument('--rmax', type=float,
                        help='Maximum radius to consider', default=10)
    parser.add_argument('--rstep', type=int,
                        help='Number of radius steps', default=200)
    args = parser.parse_args()

    r = np.linspace(args.r0, args.rmax, args.rstep+1)
    rho = woods_saxon(r, args.nucleons)
    rho0 = args.nucleons/(np.sum(r**2*rho)*4*np.pi*(r[1]-r[0]))
    rho *= rho0
    rho2 = fb_model(r, FB_6_6A)
    rho3 = fb_model(r, FB_6_6B)
    rho4 = sog_model(r, SOG_6_6)*12
    print(np.array([r, rho/2, rho2, rho4]).T)
    print(np.sum(r**2*rho)*4*np.pi*(r[1]-r[0]))
    print(np.sum(r**2*rho2)*4*np.pi*(r[1]-r[0]))
    print(np.sum(r**2*rho3)*4*np.pi*(r[1]-r[0]))
    print(np.sum(r**2*rho4)*4*np.pi*(r[1]-r[0]))


if __name__ == '__main__':
    main()
