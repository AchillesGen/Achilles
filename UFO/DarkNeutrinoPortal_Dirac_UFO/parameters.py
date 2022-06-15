# This file was automatically created by FeynRules 2.3.36
# Mathematica version: 12.0.0 for Microsoft Windows (64-bit) (April 6, 2019)
# Date: Tue 10 Aug 2021 11:00:24


from __future__ import absolute_import

from .object_library import all_parameters, Parameter


from .function_library import complexconjugate, re, im, csc, sec, acsc, asec, cot

# This is a default parameter object representing 0.
ZERO = Parameter(name = 'ZERO',
                 nature = 'internal',
                 type = 'real',
                 value = '0.0',
                 texname = '0')

# User-defined parameters.
aEWM1 = Parameter(name = 'aEWM1',
                  nature = 'external',
                  type = 'real',
                  value = 127.9,
                  texname = '\\text{aEWM1}',
                  lhablock = 'SMINPUTS',
                  lhacode = [ 1 ])

Gf = Parameter(name = 'Gf',
               nature = 'external',
               type = 'real',
               value = 0.0000116637,
               texname = 'G_f',
               lhablock = 'SMINPUTS',
               lhacode = [ 2 ])

aS = Parameter(name = 'aS',
               nature = 'external',
               type = 'real',
               value = 0.1184,
               texname = '\\alpha _s',
               lhablock = 'SMINPUTS',
               lhacode = [ 3 ])

ymdo = Parameter(name = 'ymdo',
                 nature = 'external',
                 type = 'real',
                 value = 0.00504,
                 texname = '\\text{ymdo}',
                 lhablock = 'YUKAWA',
                 lhacode = [ 1 ])

ymup = Parameter(name = 'ymup',
                 nature = 'external',
                 type = 'real',
                 value = 0.00255,
                 texname = '\\text{ymup}',
                 lhablock = 'YUKAWA',
                 lhacode = [ 2 ])

yms = Parameter(name = 'yms',
                nature = 'external',
                type = 'real',
                value = 0.101,
                texname = '\\text{yms}',
                lhablock = 'YUKAWA',
                lhacode = [ 3 ])

ymc = Parameter(name = 'ymc',
                nature = 'external',
                type = 'real',
                value = 1.27,
                texname = '\\text{ymc}',
                lhablock = 'YUKAWA',
                lhacode = [ 4 ])

ymb = Parameter(name = 'ymb',
                nature = 'external',
                type = 'real',
                value = 4.7,
                texname = '\\text{ymb}',
                lhablock = 'YUKAWA',
                lhacode = [ 5 ])

ymt = Parameter(name = 'ymt',
                nature = 'external',
                type = 'real',
                value = 172,
                texname = '\\text{ymt}',
                lhablock = 'YUKAWA',
                lhacode = [ 6 ])

yme = Parameter(name = 'yme',
                nature = 'external',
                type = 'real',
                value = 0.000511,
                texname = '\\text{yme}',
                lhablock = 'YUKAWA',
                lhacode = [ 11 ])

ymm = Parameter(name = 'ymm',
                nature = 'external',
                type = 'real',
                value = 0.10566,
                texname = '\\text{ymm}',
                lhablock = 'YUKAWA',
                lhacode = [ 13 ])

ymtau = Parameter(name = 'ymtau',
                  nature = 'external',
                  type = 'real',
                  value = 1.777,
                  texname = '\\text{ymtau}',
                  lhablock = 'YUKAWA',
                  lhacode = [ 15 ])

aD = Parameter(name = 'aD',
               nature = 'external',
               type = 'real',
               value = 0.25,
               texname = '\\text{aD}',
               lhablock = 'FRBlock',
               lhacode = [ 1 ])

aEMeps2 = Parameter(name = 'aEMeps2',
                    nature = 'external',
                    type = 'real',
                    value = 2.e-10,
                    texname = '\\text{aEMeps2}',
                    lhablock = 'FRBlock',
                    lhacode = [ 2 ])

PMNSmu42 = Parameter(name = 'PMNSmu42',
                     nature = 'external',
                     type = 'real',
                     value = 9.e-7,
                     texname = '\\text{PMNSmu42}',
                     lhablock = 'FRBlock',
                     lhacode = [ 3 ])

MZ = Parameter(name = 'MZ',
               nature = 'external',
               type = 'real',
               value = 91.1876,
               texname = '\\text{MZ}',
               lhablock = 'MASS',
               lhacode = [ 23 ])

Me = Parameter(name = 'Me',
               nature = 'external',
               type = 'real',
               value = 0.000511,
               texname = '\\text{Me}',
               lhablock = 'MASS',
               lhacode = [ 11 ])

MMU = Parameter(name = 'MMU',
                nature = 'external',
                type = 'real',
                value = 0.10566,
                texname = '\\text{MMU}',
                lhablock = 'MASS',
                lhacode = [ 13 ])

MTA = Parameter(name = 'MTA',
                nature = 'external',
                type = 'real',
                value = 1.777,
                texname = '\\text{MTA}',
                lhablock = 'MASS',
                lhacode = [ 15 ])

MU = Parameter(name = 'MU',
               nature = 'external',
               type = 'real',
               value = 0.00255,
               texname = 'M',
               lhablock = 'MASS',
               lhacode = [ 2 ])

MC = Parameter(name = 'MC',
               nature = 'external',
               type = 'real',
               value = 1.27,
               texname = '\\text{MC}',
               lhablock = 'MASS',
               lhacode = [ 4 ])

MT = Parameter(name = 'MT',
               nature = 'external',
               type = 'real',
               value = 172,
               texname = '\\text{MT}',
               lhablock = 'MASS',
               lhacode = [ 6 ])

MD = Parameter(name = 'MD',
               nature = 'external',
               type = 'real',
               value = 0.00504,
               texname = '\\text{MD}',
               lhablock = 'MASS',
               lhacode = [ 1 ])

MS = Parameter(name = 'MS',
               nature = 'external',
               type = 'real',
               value = 0.101,
               texname = '\\text{MS}',
               lhablock = 'MASS',
               lhacode = [ 3 ])

MB = Parameter(name = 'MB',
               nature = 'external',
               type = 'real',
               value = 4.7,
               texname = '\\text{MB}',
               lhablock = 'MASS',
               lhacode = [ 5 ])

MH = Parameter(name = 'MH',
               nature = 'external',
               type = 'real',
               value = 125,
               texname = '\\text{MH}',
               lhablock = 'MASS',
               lhacode = [ 25 ])

MZD = Parameter(name = 'MZD',
                nature = 'external',
                type = 'real',
                value = 0.03,
                texname = '\\text{MZD}',
                lhablock = 'MASS',
                lhacode = [ 9000005 ])

MvD = Parameter(name = 'MvD',
                nature = 'external',
                type = 'real',
                value = 0.42,
                texname = '\\text{MvD}',
                lhablock = 'MASS',
                lhacode = [ 9000006 ])

Mproton = Parameter(name = 'Mproton',
                    nature = 'external',
                    type = 'real',
                    value = 0.938,
                    texname = '\\text{M_{proton}}',
                    lhablock = 'MASS',
                    lhacode = [ 2212 ])

Mneutron = Parameter(name = 'Mneutron',
                    nature = 'external',
                    type = 'real',
                    value = 0.938,
                    texname = '\\text{M_{neutron}}',
                    lhablock = 'MASS',
                    lhacode = [ 2112 ])

Mcarbon = Parameter(name = 'Mcarbon',
                    nature = 'external',
                    type = 'real',
                    value = 11.178,
                    texname = '\\text{M_{carbon}}',
                    lhablock = 'MASS',
                    lhacode = [ 1000060120 ])

Margon = Parameter(name = 'Margon',
                   nature = 'external',
                   type = 'real',
                   value = 39.948,
                   texname = '\\text{M_{argon}}',
                   lhablock = 'MASS',
                   lhacode = [ 1000180400 ])

WZ = Parameter(name = 'WZ',
               nature = 'external',
               type = 'real',
               value = 2.4952,
               texname = '\\text{WZ}',
               lhablock = 'DECAY',
               lhacode = [ 23 ])

WW = Parameter(name = 'WW',
               nature = 'external',
               type = 'real',
               value = 2.085,
               texname = '\\text{WW}',
               lhablock = 'DECAY',
               lhacode = [ 24 ])

WT = Parameter(name = 'WT',
               nature = 'external',
               type = 'real',
               value = 1.50833649,
               texname = '\\text{WT}',
               lhablock = 'DECAY',
               lhacode = [ 6 ])

WH = Parameter(name = 'WH',
               nature = 'external',
               type = 'real',
               value = 0.00407,
               texname = '\\text{WH}',
               lhablock = 'DECAY',
               lhacode = [ 25 ])

WZD = Parameter(name = 'WZD',
                nature = 'external',
                type = 'real',
                value = 5.329174407586219e-12,
                texname = '\\text{WZD}',
                lhablock = 'DECAY',
                lhacode = [ 9000005 ])

WvD = Parameter(name = 'WvD',
                nature = 'external',
                type = 'real',
                value = 4.630139622813411e-6,
                texname = '\\text{WvD}',
                lhablock = 'DECAY',
                lhacode = [ 9000006 ])

PMNSmu4 = Parameter(name = 'PMNSmu4',
                    nature = 'internal',
                    type = 'real',
                    value = 'cmath.sqrt(PMNSmu42)',
                    texname = '\\text{PMNSmu4}')

aEW = Parameter(name = 'aEW',
                nature = 'internal',
                type = 'real',
                value = '1/aEWM1',
                texname = '\\alpha _{\\text{EW}}')

eeps = Parameter(name = 'eeps',
                 nature = 'internal',
                 type = 'real',
                 value = '2*cmath.sqrt(aEMeps2)*cmath.sqrt(cmath.pi)',
                 texname = '\\text{eeps}')

G = Parameter(name = 'G',
              nature = 'internal',
              type = 'real',
              value = '2*cmath.sqrt(aS)*cmath.sqrt(cmath.pi)',
              texname = 'G')

gD = Parameter(name = 'gD',
               nature = 'internal',
               type = 'real',
               value = '2*cmath.sqrt(aD)*cmath.sqrt(cmath.pi)',
               texname = '\\text{gD}')

MW = Parameter(name = 'MW',
               nature = 'internal',
               type = 'real',
               value = 'cmath.sqrt(MZ**2/2. + cmath.sqrt(MZ**4/4. - (aEW*cmath.pi*MZ**2)/(Gf*cmath.sqrt(2))))',
               texname = 'M_W')

ee = Parameter(name = 'ee',
               nature = 'internal',
               type = 'real',
               value = '2*cmath.sqrt(aEW)*cmath.sqrt(cmath.pi)',
               texname = 'e')

sw2 = Parameter(name = 'sw2',
                nature = 'internal',
                type = 'real',
                value = '1 - MW**2/MZ**2',
                texname = '\\text{sw2}')

cw = Parameter(name = 'cw',
               nature = 'internal',
               type = 'real',
               value = 'cmath.sqrt(1 - sw2)',
               texname = 'c_w')

sw = Parameter(name = 'sw',
               nature = 'internal',
               type = 'real',
               value = 'cmath.sqrt(sw2)',
               texname = 's_w')

g1 = Parameter(name = 'g1',
               nature = 'internal',
               type = 'real',
               value = 'ee/cw',
               texname = 'g_1')

gw = Parameter(name = 'gw',
               nature = 'internal',
               type = 'real',
               value = 'ee/sw',
               texname = 'g_w')

vev = Parameter(name = 'vev',
                nature = 'internal',
                type = 'real',
                value = '(2*MW*sw)/ee',
                texname = '\\text{vev}')

lam = Parameter(name = 'lam',
                nature = 'internal',
                type = 'real',
                value = 'MH**2/(2.*vev**2)',
                texname = '\\text{lam}')

yb = Parameter(name = 'yb',
               nature = 'internal',
               type = 'real',
               value = '(ymb*cmath.sqrt(2))/vev',
               texname = '\\text{yb}')

yc = Parameter(name = 'yc',
               nature = 'internal',
               type = 'real',
               value = '(ymc*cmath.sqrt(2))/vev',
               texname = '\\text{yc}')

ydo = Parameter(name = 'ydo',
                nature = 'internal',
                type = 'real',
                value = '(ymdo*cmath.sqrt(2))/vev',
                texname = '\\text{ydo}')

ye = Parameter(name = 'ye',
               nature = 'internal',
               type = 'real',
               value = '(yme*cmath.sqrt(2))/vev',
               texname = '\\text{ye}')

ym = Parameter(name = 'ym',
               nature = 'internal',
               type = 'real',
               value = '(ymm*cmath.sqrt(2))/vev',
               texname = '\\text{ym}')

ys = Parameter(name = 'ys',
               nature = 'internal',
               type = 'real',
               value = '(yms*cmath.sqrt(2))/vev',
               texname = '\\text{ys}')

yt = Parameter(name = 'yt',
               nature = 'internal',
               type = 'real',
               value = '(ymt*cmath.sqrt(2))/vev',
               texname = '\\text{yt}')

ytau = Parameter(name = 'ytau',
                 nature = 'internal',
                 type = 'real',
                 value = '(ymtau*cmath.sqrt(2))/vev',
                 texname = '\\text{ytau}')

yup = Parameter(name = 'yup',
                nature = 'internal',
                type = 'real',
                value = '(ymup*cmath.sqrt(2))/vev',
                texname = '\\text{yup}')

muH = Parameter(name = 'muH',
                nature = 'internal',
                type = 'real',
                value = 'cmath.sqrt(lam*vev**2)',
                texname = '\\mu')

I1b11 = Parameter(name = 'I1b11',
                  nature = 'internal',
                  type = 'complex',
                  value = 'ydo',
                  texname = '\\text{I1b11}')

I1b22 = Parameter(name = 'I1b22',
                  nature = 'internal',
                  type = 'complex',
                  value = 'ys',
                  texname = '\\text{I1b22}')

I1b33 = Parameter(name = 'I1b33',
                  nature = 'internal',
                  type = 'complex',
                  value = 'yb',
                  texname = '\\text{I1b33}')

I2b11 = Parameter(name = 'I2b11',
                  nature = 'internal',
                  type = 'complex',
                  value = 'yup',
                  texname = '\\text{I2b11}')

I2b22 = Parameter(name = 'I2b22',
                  nature = 'internal',
                  type = 'complex',
                  value = 'yc',
                  texname = '\\text{I2b22}')

I2b33 = Parameter(name = 'I2b33',
                  nature = 'internal',
                  type = 'complex',
                  value = 'yt',
                  texname = '\\text{I2b33}')

I3b11 = Parameter(name = 'I3b11',
                  nature = 'internal',
                  type = 'complex',
                  value = 'yup',
                  texname = '\\text{I3b11}')

I3b22 = Parameter(name = 'I3b22',
                  nature = 'internal',
                  type = 'complex',
                  value = 'yc',
                  texname = '\\text{I3b22}')

I3b33 = Parameter(name = 'I3b33',
                  nature = 'internal',
                  type = 'complex',
                  value = 'yt',
                  texname = '\\text{I3b33}')

I4b11 = Parameter(name = 'I4b11',
                  nature = 'internal',
                  type = 'complex',
                  value = 'ydo',
                  texname = '\\text{I4b11}')

I4b22 = Parameter(name = 'I4b22',
                  nature = 'internal',
                  type = 'complex',
                  value = 'ys',
                  texname = '\\text{I4b22}')

I4b33 = Parameter(name = 'I4b33',
                  nature = 'internal',
                  type = 'complex',
                  value = 'yb',
                  texname = '\\text{I4b33}')

