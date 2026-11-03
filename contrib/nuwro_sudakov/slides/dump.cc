// dump NuWro nucleus inputs for 12C: r, rho [fm^-3], kF_p, kF_n [MeV], Ef [MeV], corr tables
#include "kaskada7.h"
#include "input_data.h"
#include <cstdio>
int main(int argc,char**argv){
  set_dirs(argv[0]);
  params p; p.read(argv[1]); p.nucleus_p=6; p.nucleus_n=6;
  nucleus *n=make_nucleus(p);
  printf("# R=%g fm\n", n->radius()/fermi);
  for(int i=0;i<=900;i++){ double r=i*0.01*fermi;
    particle pr(PDG::pdg_proton,PDG::mass_proton); pr.r=vec(r,0,0);
    printf("%g %g %g %g\n", r/fermi, n->density(r)*fermi3, n->localkf(pr), n->Ef(pr)); }
  input_data in; in.initialize(p); in.load_data();
  data_container* c=in.get_nucl_data_container(0,6,6);
  FILE*f=fopen("corr.txt","w");
  for(int i=0;i<=1000;i++){ double s=i*0.01*fermi; c->set_input_point(s);
    fprintf(f,"%g %g %g %g %g\n", s/fermi, c->get_value(1), c->get_value(2), c->get_value(3), c->get_value(4)); }
  fclose(f);
}
