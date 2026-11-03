// Dump per-event summary from a NuWro kaskada ROOT file.
// columns: Ek_in[MeV] noi nod0 nod2 nod3 nod4..8sum firstx firsty firstz startx starty startz npost
#include "event1.h"
#include <TFile.h>
#include <TTree.h>
#include <cstdio>
int main(int argc,char**argv){
  TFile f(argv[1]); TTree*t=(TTree*)f.Get("treeout"); event*e=new event; t->SetBranchAddress("e",&e);
  long n=t->GetEntries();
  for(long i=0;i<n;i++){ t->GetEntry(i);
    int noi=e->number_of_interactions();
    double fx=0,fy=0,fz=0; 
    if(noi>0 && e->all.size()>0){ fx=e->all[0].r.x; fy=e->all[0].r.y; fz=e->all[0].r.z;}
    int pions=0; for(int k=4;k<=8;k++) pions+=e->nod[k];
    printf("%.4f %d %d %d %d %d %.5f %.5f %.5f %.5f %.5f %.5f %d\n", e->out[0].Ek(), noi, e->nod[0], e->nod[2], e->nod[3], pions,
      fx,fy,fz, e->out[0].r.x,e->out[0].r.y,e->out[0].r.z,(int)e->post.size());
  }
}
