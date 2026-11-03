import sys, json, os, csv; sys.path.insert(0,'/home/claude/work/ana')
from common import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':15,'axes.grid':True,'grid.alpha':0.3,'axes.spines.top':False,'axes.spines.right':False,'legend.frameon':False,'mathtext.fontset':'cm'})
OUT='/home/claude/slides/figs'; CA,CB,CX='#1f77b4','#d62728','#2ca02c'
def save(f,n): f.tight_layout(); f.savefig(f'{OUT}/{n}.pdf'); plt.close(f)
def readcsv(f):
    rows=[r for r in csv.reader(l for l in open(f) if not l.startswith('#'))]; return np.array([[float(x) for x in r] for r in rows[1:]])
dR=readcsv('/home/claude/work/data/pC_reaction_xsec_Carlson.csv'); dT=readcsv('/home/claude/work/data/C12_transparency.csv')
s2=json.load(open('/home/claude/work/ana/step2.json')); sb=load('stepb',''); chi=load('chi','chi_'); noi=json.load(open('/home/claude/work/ana/noi.json'))
# step dependence 2x2
fig,axs=plt.subplots(2,2,figsize=(10,7.2))
for ax,(obs,E) in zip(axs.flat,[('sR',10),('sR',300),('T',300),('T',1000)]):
    r=s2[f'{obs}_{E}']; x=np.array(r['steps']); y=np.array(r['vals']); e=np.array(r['errs'])
    ax.errorbar(x,y,e,fmt='o',color=CA,ms=5,label='(a) at step $\\Delta$')
    m=x<=0.8; c,_,_,_=wfit(x[m],y[m],e[m],2); xx=np.linspace(0,0.8,100); ax.plot(xx,c[0]+c[1]*xx+c[2]*xx**2,color=CA,lw=1)
    ax.errorbar([0],[r['zero_q']],[r['zero_q_e']],fmt='s',color=CX,ms=7,label='(a) extrapolated to $\\Delta\\to0$')
    v,ev=(sigR if obs=='sR' else trans)(sb[f'{obs}_b_g2.0_E{E}']); ax.errorbar([-0.05],[v],[ev],fmt='D',color=CB,ms=7,label='(b), no step')
    ax.axvline(0.2,color='gray',ls=':'); ax.text(0.22,ax.get_ylim()[0] if False else y.min(),'default',color='gray',fontsize=10)
    ax.set_title(('$\\sigma_R$ [mb], ' if obs=='sR' else 'transparency, ')+f'{E} MeV',fontsize=13); ax.set_xlabel('$\\Delta$ [fm]')
axs[0,0].legend(fontsize=10)
save(fig,'res_step')
# data comparison
fig,axs=plt.subplots(1,2,figsize=(11,4.2))
for ax,(obs,fn,dd,yl) in zip(axs,[('sR',sigR,dR,'$\\sigma_R$ [mb]'),('T',trans,dT,'transparency')]):
    ax.errorbar(dd[:,0],dd[:,1],dd[:,2],fmt='o',color='k',ms=3,elinewidth=0.8,label='data (digitised)')
    for v,lab,ls in [('a','(a), $\\Delta=0.2$ fm','-'),('b','(b)','--')]:
        p=np.array(sorted((s['Ek'],)+fn(s) for t,s in chi.items() if t.startswith(f'chi_{obs}_{v}_')))
        ax.plot(p[:,0],p[:,1],ls,color=CA if v=='a' else CB,lw=2,label=lab)
    ax.set_xscale('log'); ax.set_ylabel(yl); ax.set_xlabel('kinetic energy [MeV]'); ax.legend(fontsize=10)
save(fig,'res_data')
fig,axs=plt.subplots(1,2,figsize=(11,3.6))
for ax,(obs,dd) in zip(axs,[('sR',dR),('T',dT)]):
    Es=[];b=[];be=[]
    for k,r in s2.items():
        if r['obs']==obs: Es.append(r['E']); b.append(r['bias02']/r['zero_q']*100); be.append(r['bias02_e']/r['zero_q']*100)
    o=np.argsort(Es)
    o2=np.argsort(dd[:,0]); ax.fill_between(dd[o2,0],-(dd[o2,2]/dd[o2,1]*100),(dd[o2,2]/dd[o2,1]*100),color='gray',alpha=0.3,step='mid',label='data relative error')
    ax.errorbar(np.array(Es)[o],np.array(b)[o],np.array(be)[o],fmt='o',color=CA,ms=5,label='(a): bias at $\\Delta=0.2$ fm')
    ax.axhline(0,color='k',lw=0.7); ax.set_xscale('log'); ax.set_ylim(-12,12); ax.set_ylabel('% of $\\Delta\\to0$ value'); ax.set_xlabel('kinetic energy [MeV]')
    ax.set_title('$\\sigma_R$' if obs=='sR' else 'transparency',fontsize=13); ax.legend(fontsize=10,loc='lower right')
save(fig,'res_bias')
# N_coll
fig,axs=plt.subplots(1,3,figsize=(12.5,3.9))
for ax,key,title in [(axs[0],'T_300','transparency mode, 300 MeV'),(axs[1],'T_1000','transparency mode, 1 GeV')]:
    r=noi[key]; ax.errorbar(r['steps'],r['vals'],r['errs'],fmt='o',ms=5,color=CA,label='(a)')
    ax.errorbar([0],[r['zero']],[r['zero_e']],fmt='s',ms=7,color=CX,label='(a), $\\Delta\\to0$'); ax.errorbar([-0.05],[r['b']],[r['b_e']],fmt='D',ms=7,color=CB,label='(b)')
    ax.axvline(0.2,color='gray',ls=':'); ax.set_title(title,fontsize=12); ax.set_xlabel('$\\Delta$ [fm]'); ax.set_ylabel('$\\langle N_{\\rm coll}\\rangle$')
axs[0].legend(fontsize=10)
ax=axs[2]; xs=[];ys=[];es=[]
for st,f in [(0.0125,'sfstep_0.0125_300.txt'),(0.025,'sfstep_0.025_300.txt'),(0.05,'sfstep_0.05_300.txt'),(0.1,'sfstep_0.1_300.txt'),(0.2,'sf_0_300.txt'),(0.4,'sfstep_0.4_300.txt'),(0.8,'sfstep_0.8_300.txt')]:
    a=np.loadtxt('/home/claude/work/sf/'+f); m=a[:,0]==0; v=a[m,2]; xs.append(st); ys.append(v.mean()); es.append(v.std()/np.sqrt(m.sum()))
ax.errorbar(xs,ys,es,fmt='o',ms=5,color=CA,label='(a)')
for i,f in enumerate(['sf_1_300.txt','sfseg_0.5_300.txt','sfseg_8_300.txt']):
    a=np.loadtxt('/home/claude/work/sf/'+f); m=a[:,0]==0; ax.errorbar([-0.03-0.03*i],[a[m,2].mean()],[a[m,2].std()/np.sqrt(m.sum())],fmt='D',ms=7,color=CB,label='(b), $L_{\\rm seg}$=2, 0.5, 8 fm' if i==0 else None)
ax.axvline(0.2,color='gray',ls=':'); ax.set_xlabel('$\\Delta$ [fm]'); ax.set_ylabel('$\\langle N_{\\rm coll}\\rangle$'); ax.set_title('QE-SF branch, $T_{N1}$=300 MeV',fontsize=12); ax.legend(fontsize=10)
save(fig,'res_ncoll')
# reweight
fig,axs=plt.subplots(1,2,figsize=(11,4))
ax=axs[0]
for cfg,col in [('T300','#9467bd'),('T1000','#8c564b')]:
    a=np.load(f'/home/claude/work/rw/rw_b_{cfg}_nom.npy').astype(float); W=np.loadtxt(f'/home/claude/work/rw/rw_b_{cfg}_nom.w.txt')[:,2:]; P0=(a[:,1]==0)*1.
    T=[(W[:,i]*P0).sum()/W[:,i].sum() for i in range(4)]
    ax.plot([0.6,0.7,0.76,0.82,1.0],[T[0],T[1],P0.mean(),T[2],T[3]],'-o',ms=6,mfc='none',color=col,label=f'reweighted from 0.76, {cfg[1:]} MeV')
    for v in [0.6,1.0]:
        b=np.load(f'/home/claude/work/rw/rw_b_{cfg}_mfp{v}.npy').astype(float); p=(b[:,1]==0).mean(); ax.errorbar([v],[p],[np.sqrt(p*(1-p)/len(b))],fmt='D',ms=8,color=col)
ax.plot([],[],'D',color='gray',label='regenerated'); ax.set_xlabel('kaskada_NN_mfp_scale'); ax.set_ylabel('transparency $P_0$'); ax.legend(fontsize=10)
ax=axs[1]; rw=json.load(open('/home/claude/work/ana/rw.json'))
VARS=["mfp=0.6","mfp=0.7","mfp=0.82","mfp=1.0","xs=0.8","xs=0.9","xs=1.1","xs=1.25"]
for alg,mk,col in [('a','o',CA),('b','D',CB)]:
    for cfg,ls in [('T300','-'),('T1000','--'),('sR50',':'),('sR300','-.')]:
        ax.plot(range(8),[rw[f'{alg}_{cfg}_{v}']['ess_frac'] for v in VARS],mk+ls,color=col,ms=5,lw=1,label=f'({alg}) {cfg}')
ax.set_xticks(range(8)); ax.set_xticklabels(VARS,rotation=45,fontsize=10); ax.set_ylabel('ESS / N'); ax.legend(fontsize=8,ncol=2)
save(fig,'res_rw')
print('ok')
