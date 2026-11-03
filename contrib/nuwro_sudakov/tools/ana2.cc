#include "event1.h"
#include <TFile.h>
#include <TTree.h>
#include <cstdio>
int main(int argc,char**argv){
  TFile f(argv[1]); TTree*t=(TTree*)f.Get("treeout"); event*e=new event; t->SetBranchAddress("e",&e);
  for(long i=0;i<t->GetEntries();i++){ t->GetEntry(i); printf("%ld %d %d %d",i,e->pr,e->nr,(int)e->post.size());
    for(auto&q:e->post) printf(" %d %.6f %.6f %.6f",q.pdg,q.x,q.y,q.z); for(int k=0;k<18;k++) if(k!=12&&k!=13&&k!=14) printf(" %d",e->nod[k]); printf("\n"); }
}
