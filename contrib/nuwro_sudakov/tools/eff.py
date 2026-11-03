import numpy as np, re
txt=open('/home/claude/src/nuwro/data/input/kaskada_NN_xsec_2.dat').read()
rows=[]; cur={}
for line in txt.splitlines():
    line=line.strip()
    if line=='-':
        if cur: rows.append(cur); cur={}
        continue
    m=re.match(r'(\w+):\s*([-\deE.+]+)',line)
    if m: cur[m.group(1)]=float(m.group(2))
if cur: rows.append(cur)
E=np.array([r['energy'] for r in rows]); ii=np.array([r.get('xsec_ii',np.nan) for r in rows]); ij=np.array([r.get('xsec_ij',np.nan) for r in rows])
def fill(y):
    m=~np.isnan(y); return np.interp(E,E[m],y[m])
ii=fill(ii); ij=fill(ij)
M=938.919; rng=np.random.default_rng(0)
print("Ek  kF  sup(sig_ij)  mean(sig_ij)  ratio   | ii")
for Ek in [30,50,100,200,500]:
  for kf in [200,263]:
    p1=np.sqrt(Ek**2+2*M*Ek); E1=Ek+M
    n=400000; q=kf*rng.uniform(0,1,n)**(1/3); c=rng.uniform(-1,1,n)
    E2=np.sqrt(q*q+M*M); Eke=(E1*E2-p1*q*c)/M-M
    s=np.interp(Eke,E,ij); s2=np.interp(Eke,E,ii)
    # sup over ball, as in the patch's bound: Ekeff range
    qq=kf; E2m=np.sqrt(qq*qq+M*M)
    Emax=(E1*E2m+p1*qq)/M-M; Emin=0 if p1<=qq else (E1*E2m-p1*qq)/M-M
    grid=np.r_[Emin,Emax,E[(E>Emin)&(E<Emax)]]
    sup=np.interp(grid,E,ij).max(); sup2=np.interp(grid,E,ii).max()
    print(f"{Ek:4d} {kf:4d} {sup:9.1f} {s.mean():9.1f} {s.mean()/sup:7.3f} | {sup2:8.1f} {s2.mean():8.1f} {s2.mean()/sup2:6.3f}")
