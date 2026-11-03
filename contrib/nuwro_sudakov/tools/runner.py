#!/usr/bin/env python3
"""Run NuWro kaskada (a) or (b) and summarise. Usage: runner.py jobs.json outdir nproc"""
import json, subprocess, sys, os, numpy as np, time
from concurrent.futures import ThreadPoolExecutor
RR="/home/claude/rr.sh"
MN=938.272  # proton mass used for beam_energy (total energy) [NuWro pdg.h mass_proton=938.27203]
def run(job, outdir):
    tag=job['tag']; out=os.path.join(outdir,tag)
    if os.path.exists(out+'.sum.json'): return json.load(open(out+'.sum.json'))
    m = 938.27203 if job.get('pdg',2212)==2212 else 939.56536
    args=[RR, job['bin'], '-i', job.get('input','/home/claude/src/nuwro_a/data/kaskada.txt'), '-o', out+'.root',
          '-p', f"number_of_events={job['n']}", '-p', f"beam_particle={job.get('pdg',2212)}",
          '-p', f"beam_energy={m+job['Ek']:.5f}", '-p', f"random_seed={job['seed']}",
          '-p', f"beam_placement={job['placement']}", '-p', f"step={job['step']}"]
    for k,v in job.get('extra',{}).items(): args += ['-p', f"{k}={v}"]
    t0=time.time(); 
    with open(out+'.log','w') as lf: subprocess.run(args, stdout=lf, stderr=subprocess.STDOUT, cwd=outdir, check=True)
    cpu=time.time()-t0
    if os.path.exists(out+'.root.sud.txt'): job=dict(job, sud=open(out+'.root.sud.txt').read())
    if os.path.exists(out+'.root.w.txt'):
        import shutil; shutil.move(out+'.root.w.txt', out+'.w.txt')
    anadir='/home/claude/tools/'+job.get('ana','a')
    txt=subprocess.run([RR, anadir+'/ana', out+'.root'], capture_output=True, text=True, cwd=anadir).stdout
    if not txt.strip(): raise RuntimeError('ana failed '+tag)
    a=np.loadtxt(txt.splitlines(), ndmin=2)
    np.save(out+'.npy', a.astype(np.float32))
    os.remove(out+'.root')
    n=len(a); k=int((a[:,1]>0).sum())
    s=dict(job, nev=n, nint=k, pint=k/n, wall=cpu)
    json.dump(s, open(out+'.sum.json','w')); return s
if __name__=='__main__':
    jobs=json.load(open(sys.argv[1])); outdir=sys.argv[2]; os.makedirs(outdir,exist_ok=True)
    with ThreadPoolExecutor(int(sys.argv[3])) as ex:
        for s in ex.map(lambda j: run(j,os.path.join('/home/claude/work',j['dir']) if 'dir' in j else outdir), jobs):
            print(s['tag'], s['nev'], s['pint'], round(s['wall'],1), flush=True)
