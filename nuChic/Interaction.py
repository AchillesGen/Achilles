from nuChic.Constants import mqe as mN, MeV, mb, GeV
import numpy as np

# HZETRN parameters for pp, pn, and nn interactions
a = 5.0 * MeV              
b = 0.199/np.sqrt(MeV)     
c = 0.451 * MeV**-0.258      
d = 25.0 * MeV               
e = 134.0 * MeV              
f = 1.187 * MeV**-0.35       
g = 0.1 * MeV                
h = 0.282 * MeV              

# PDG parameters for pp, pn, and nn interactions
Zpp  = 33.45 * mb            
Zpn  = 35.80 * mb              
Y1pp = 42.53 * mb              
Y1pn = 40.15 * mb              
Y2pp = 33.34 * mb              
Y2pn = 30.00 * mb              
B    = 0.308 * mb              
s1   = 1.0   * GeV**2           
s0   = (5.38 * GeV)**2 
n1   = 0.458
n2   = 0.545

# JWN parameters for pp, pn, and nn interactions
gamma = 52.5 * mb * GeV**(0.16)
alpha = 0.00369 / MeV
beta  = 0.00895741 * MeV**(-0.8)

def sigma_pp(plab):
    Tlab = np.sqrt(plab**2+mN**2)-mN
    if(plab < 1.8):
        if(Tlab >= 0.025):
            return (1+a/Tlab)*(40+109*np.cos(b*np.sqrt(Tlab))*np.exp(-c*(Tlab-d)**(0.258)))
        else:
            return np.exp(6.51*np.exp(-(Tlab/e)**(0.7)))
    elif(plab <= 4.7):
        return gamma/plab**0.16
    else:
        s = 2*mN*(mN+np.sqrt(plab**2+mN**2))
        return Zpp + B*np.log(s/s0)**2 + Y1pp*(s1/s)**n1 - Y2pp*(s1/s)**n2
    
def sigma_np(plab):
    Tlab = np.sqrt(plab**2+mN**2)-mN
    if(plab < 0.5):
        if(Tlab >= 0.1 * MeV):
            return 38 + 12500*np.exp(-f*(Tlab-g)**0.35)
        else:
            return 26000 * np.exp(-(Tlab/h)**0.3)
    elif(plab <= 2.0):
        return 40 + 10*np.cos(alpha*plab - 0.943)*np.exp(-beta*plab**0.8+2)
    else:
        s = 2*mN*(mN+np.sqrt(plab**2+mN**2))
        return Zpn + B*np.log(s/s0)**2 + Y1pn*(s1/s)**n1 - Y2pn*(s1/s)**n2
