#include "event1.h"
#include <TFile.h>
#include <TTree.h>
#include <cstdio>
int main(int argc,char**argv){
  TFile f(argv[1]); TTree*t=(TTree*)f.Get("treeout"); event*e=new event; t->SetBranchAddress("e",&e);
  long n=t->GetEntries(), q=0, tr=0, co=0;
  for(long i=0;i<n;i++){ t->GetEntry(i); q+=e->flag.qel!=0; tr+=e->flag.isTransparent!=0; co+=e->flag.isCorrelated!=0;}
  printf("%ld events: qel=%ld transparent=%ld correlated=%ld\n",n,q,tr,co);
}
