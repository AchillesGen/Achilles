import sys; sys.path.insert(0,'/home/claude/work/ana')
from common import *
import csv, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':9,'axes.grid':True,'grid.alpha':0.3,'figure.dpi':150})
OUT='/home/claude/work/report/figs'; os.makedirs(OUT,exist_ok=True)
C={'a':'#1f77b4','b':'#d62728','data':'k','ext':'#2ca02c'}
def readcsv(f):
    rows=[r for r in csv.reader(l for l in open(f) if not l.startswith('#'))]
    return np.array([[float(x) for x in r] for r in rows[1:]])
dR=readcsv('/home/claude/work/data/pC_reaction_xsec_Carlson.csv'); dT=readcsv('/home/claude/work/data/C12_transparency.csv')
s2=json.load(open('/home/claude/work/ana/step2.json'))
sb=load('stepb',''); chi=load('chi','chi_')

# ---- Fig 1: step dependence at selected energies, with (b) at zero step
sel=[('sR',10),('sR',30),('sR',100),('sR',300),('T',300),('T',1000)]
fig,axs=plt.subplots(2,3,figsize=(9,5.2))
for ax,(obs,E) in zip(axs.flat,sel):
    r=s2.get(f'{obs}_{E}')
    if not r: continue
    x=np.array(r['steps']); y=np.array(r['vals']); e=np.array(r['errs'])
    ax.errorbar(x,y,e,fmt='o',color=C['a'],ms=3,label='(a) runs')
    xx=np.linspace(0,max(x),100); c=[r['zero_q'],0,0]
    # refit quad for drawing
    mq=x<=0.8; cq,_,_,_=wfit(x[mq],y[mq],e[mq],2); xx=np.linspace(0,0.8,100); ax.plot(xx,cq[0]+cq[1]*xx+cq[2]*xx**2,color=C['a'],lw=0.8)
    ax.errorbar([0],[r['zero_q']],[r['zero_q_e']],fmt='s',color=C['ext'],ms=4,label='(a) extrapolated')
    t=f"{obs}_b_g2.0_E{E}"
    if t in sb:
        v,ev=(sigR if obs=='sR' else trans)(sb[t]); ax.errorbar([-0.03],[v],[ev],fmt='D',color=C['b'],ms=4,label='(b), $L_{seg}$=2 fm')
    ax.axvline(0.2,color='gray',ls=':',lw=0.8)
    ax.set_title(f"{'$\\sigma_R$ [mb]' if obs=='sR' else 'transparency'}, $T$={E} MeV",fontsize=9)
    ax.set_xlabel('max_step $\\Delta$ [fm]')
axs[0,0].legend(fontsize=7)
fig.tight_layout(); fig.savefig(f'{OUT}/fig1_step.png'); plt.close(fig)

# ---- Fig 2: sigma_R and T vs energy: data, (a) default, (b); bias panel vs data errors
fig,axs=plt.subplots(2,2,figsize=(9,6),gridspec_kw={'height_ratios':[2,1]},sharex='col')
for col,(obs,fn,dd) in enumerate([('sR',sigR,dR),('T',trans,dT)]):
    ax=axs[0,col]
    ax.errorbar(dd[:,0],dd[:,1],dd[:,2],fmt='o',color='k',ms=2.5,elinewidth=0.7,label='data (digitised)')
    for v,lab in [('a','(a), $\\Delta$=0.2 fm'),('b','(b)')]:
        pts=sorted((s['Ek'],)+fn(s) for t,s in chi.items() if t.startswith(f'chi_{obs}_{v}_'))
        if pts:
            p=np.array(pts); ax.plot(p[:,0],p[:,1],color=C[v],lw=1.2 if v=='a' else 0.9,ls='-' if v=='a' else '--',label=lab)
    ax.set_xscale('log'); ax.set_ylabel('$\\sigma_R$ [mb]' if obs=='sR' else 'transparency'); ax.legend(fontsize=7)
    ax=axs[1,col]
    Es=[];b=[];be=[]
    for k,r in s2.items():
        if r['obs']==obs: Es.append(r['E']); b.append(r['bias02']/r['zero_q']*100); be.append(r['bias02_e']/r['zero_q']*100)
    o=np.argsort(Es); Es=np.array(Es)[o]; b=np.array(b)[o]; be=np.array(be)[o]
    ax.errorbar(Es,b,be,fmt='o',color=C['a'],ms=3,label='(a) bias at $\\Delta$=0.2 fm')
    ax.fill_between(dd[:,0][np.argsort(dd[:,0])],-(dd[:,2]/dd[:,1]*100)[np.argsort(dd[:,0])],(dd[:,2]/dd[:,1]*100)[np.argsort(dd[:,0])],color='gray',alpha=0.25,step='mid',label='data rel. error')
    ax.axhline(0,color='k',lw=0.6); ax.set_ylabel('% of zero-step value'); ax.set_xlabel('kinetic energy [MeV]'); ax.legend(fontsize=7)
    ax.set_ylim(-15,15)
fig.tight_layout(); fig.savefig(f'{OUT}/fig2_data.png'); plt.close(fig)
print('figs written')

# ---- Fig 3: <N_int> vs step (cascade-level bias) and SF-mode test
noi=json.load(open('/home/claude/work/ana/noi.json'))
fig,axs=plt.subplots(1,3,figsize=(9,3))
for ax,key,title in [(axs[0],'T_300','transparency mode, $T_p$=300 MeV'),(axs[1],'T_1000','transparency mode, $T_p$=1 GeV')]:
    r=noi[key]; x=np.array(r['steps']); y=np.array(r['vals']); e=np.array(r['errs'])
    ax.errorbar(x,y,e,fmt='o',ms=3,color=C['a'],label='(a)')
    ax.errorbar([0],[r['zero']],[r['zero_e']],fmt='s',ms=4,color=C['ext'],label='(a) extrapolated')
    ax.errorbar([-0.03],[r['b']],[r['b_e']],fmt='D',ms=4,color=C['b'],label='(b)')
    ax.axvline(0.2,color='gray',ls=':'); ax.set_title(title,fontsize=8); ax.set_xlabel('$\\Delta$ [fm]'); ax.set_ylabel('$\\langle N_{\\rm coll}\\rangle$')
axs[0].legend(fontsize=7)
ax=axs[2]
sf=[(0.0125,'sfstep_0.0125_300.txt'),(0.025,'sfstep_0.025_300.txt'),(0.05,'sfstep_0.05_300.txt'),(0.1,'sfstep_0.1_300.txt'),(0.2,'sf_0_300.txt'),(0.4,'sfstep_0.4_300.txt'),(0.8,'sfstep_0.8_300.txt')]
xs=[];ys=[];es=[]
for st,f in sf:
    a=np.loadtxt('/home/claude/work/sf/'+f); m=a[:,0]==0; v=a[m,2]; xs.append(st); ys.append(v.mean()); es.append(v.std()/np.sqrt(m.sum()))
ax.errorbar(xs,ys,es,fmt='o',ms=3,color=C['a'],label='(a)')
bb=[]
for f in ['sf_1_300.txt','sfseg_0.5_300.txt','sfseg_8_300.txt']:
    p='/home/claude/work/sf/'+f
    if os.path.exists(p):
        a=np.loadtxt(p); m=a[:,0]==0; bb.append((a[m,2].mean(),a[m,2].std()/np.sqrt(m.sum())))
for i,(v,e) in enumerate(bb): ax.errorbar([-0.02-0.015*i],[v],[e],fmt='D',ms=4,color=C['b'],label='(b), $L_{seg}$=2,0.5,8 fm' if i==0 else None)
ax.axvline(0.2,color='gray',ls=':'); ax.set_xlabel('$\\Delta$ [fm]'); ax.set_ylabel('$\\langle N_{\\rm coll}\\rangle$'); ax.set_title('QE-SF branch, $T_{N1}$=300 MeV',fontsize=8); ax.legend(fontsize=7)
fig.tight_layout(); fig.savefig(f'{OUT}/fig3_ncoll.png'); plt.close(fig)

# ---- Fig 4: reweighting closure
rw=json.load(open('/home/claude/work/ana/rw.json'))
fig,axs=plt.subplots(1,3,figsize=(10,3.2))
# panel 1: transparency vs mfp scale from weights (b), regenerated points
ax=axs[0]
for cfg,col in [('T300','#9467bd'),('T1000','#8c564b')]:
    xs=[0.76];ys=[];
    nom=None
    pts=[]
    for k,r in rw.items():
        if r['alg']=='b' and r['cfg']==cfg and r['var'].startswith('mfp'):
            pass
    a=np.load(f'/home/claude/work/rw/rw_b_{cfg}_nom.npy').astype(float); W=np.loadtxt(f'/home/claude/work/rw/rw_b_{cfg}_nom.w.txt')[:,2:]
    P0=(a[:,1]==0)*1.
    vals=[0.6,0.7,0.82,1.0]; T=[(W[:,i]*P0).sum()/W[:,i].sum() for i in range(4)]
    xx=[0.6,0.7,0.76,0.82,1.0]; yy=[T[0],T[1],P0.mean(),T[2],T[3]]
    ax.plot(xx,yy,'-o',ms=3,mfc='none',color=col,label=f'(b) reweighted, $T_p$={cfg[1:]}')
    for v in [0.6,1.0]:
        b=np.load(f'/home/claude/work/rw/rw_b_{cfg}_mfp{v}.npy').astype(float); p=(b[:,1]==0).mean(); e=np.sqrt(p*(1-p)/len(b))
        ax.errorbar([v],[p],[e],fmt='D',ms=4,color=col)
ax.set_xlabel('kaskada_NN_mfp_scale'); ax.set_ylabel('transparency $P_0$'); ax.legend(fontsize=7); ax.set_title('filled: regenerated',fontsize=8)
# panel 2: pulls (reweighted - regenerated)/sigma for all observables
ax=axs[1]; i=0; labs=[]
for alg,mk,col in [('b','D',C['b']),('a','o',C['a'])]:
    for k,r in sorted(rw.items()):
        if r['alg']!=alg or '<N_int>' not in r: continue
        for j,o in enumerate(['P0 (no interaction)','<N_int>','P(N_int>=2)','<N_out>','<N_elastic>']):
            ax.plot(r[o]['z'], i, mk, ms=3, color=col, mfc='none' if j else col)
        labs.append(f"({alg}) {r['cfg']} {r['var']}"); i+=1
ax.set_yticks(range(len(labs))); ax.set_yticklabels(labs,fontsize=5); ax.axvline(0,color='k',lw=0.6)
for z in (-2,2): ax.axvline(z,color='gray',ls=':',lw=0.6)
ax.set_xlabel('(reweighted $-$ regenerated)/$\\sigma$'); ax.set_title('5 observables per row',fontsize=8)
# panel 3: ESS
ax=axs[2]
VARS=["mfp=0.6","mfp=0.7","mfp=0.82","mfp=1.0","xs=0.8","xs=0.9","xs=1.1","xs=1.25"]
for alg,mk,col in [('b','D',C['b']),('a','o',C['a'])]:
    for cfg,ls in [('T300','-'),('T1000','--'),('sR50',':'),('sR300','-.')]:
        e=[rw.get(f'{alg}_{cfg}_{v}',{}).get('ess_frac',np.nan) for v in VARS]
        ax.plot(range(len(VARS)),e,mk+ls,ms=3,color=col,lw=0.8,label=f'({alg}) {cfg}')
ax.set_xticks(range(len(VARS))); ax.set_xticklabels(VARS,rotation=60,fontsize=6); ax.set_ylabel('ESS / N'); ax.legend(fontsize=5,ncol=2)
fig.tight_layout(); fig.savefig(f'{OUT}/fig4_reweight.png'); plt.close(fig)
print('fig 4 written')
