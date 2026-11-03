// Exercise the QE spectral-function branches of kaskada::kaskadaevent (transparent/correlated,
// with N1 redraws) for (a) and (b).  Output per event: flags, noi, N1-interacted proxy, npost protons, nod[0..3], nod[8].
#include "kaskada7.h"
#include "input_data.h"
#include "args.h"
#include <cstdio>
int main(int argc, char** argv)
{
  // argv: params_file sudakov nev seed T_N1[MeV]
  set_dirs(argv[0]);
  params p; p.read(argv[1]);
  p.kaskada_sudakov = atoi(argv[2]);
  int nev = atoi(argv[3]);
  p.random_seed = atoi(argv[4]);
  double T = atof(argv[5]);
  if (argc > 6) p.step = atof(argv[6]);
  if (argc > 7) p.kaskada_sudakov_segment = atof(argv[7]);
  p.sf_method = 1; p.nucleus_target = 2; p.beam_placement = 0;
  p.nucleus_p = 6; p.nucleus_n = 6;
  frandom_init(p.random_seed);
  input_data input; input.initialize(p); input.load_data();
  nucleus *nucl = make_nucleus(p);
  for (int i = 0; i < nev; i++)
  {
    event *e = new event; e->par = p; e->weight = 1;
    e->flag.qel = 1; e->flag.isTransparent = frandom() < 0.3; e->flag.isCorrelated = frandom() < 0.4;
    vec r = nucl->get_random_r() * rand_dir();
    particle nu(PDG::pdg_nu_mu, 0); nu.set_momentum(vec(0,0,1000)); nu.r = r;
    particle mu(PDG::pdg_mu, PDG::mass_mu); mu.set_momentum(vec(0,0,500)); mu.r = r;
    particle N1(PDG::pdg_proton, PDG::mass_proton);
    double pp = sqrt(T*T + 2*T*PDG::mass_proton);
    vec d = rand_dir(); N1.set_momentum(pp*d); N1.r = r;
    particle Nt(PDG::pdg_neutron, PDG::mass_neutron); Nt.r = r;
    e->in.push_back(nu); e->in.push_back(Nt);
    e->out.push_back(mu); e->out.push_back(N1);
    if (e->flag.isCorrelated)
    {
      particle N2(PDG::pdg_neutron, PDG::mass_neutron); N2.set_momentum(-300*d + 50*rand_dir()); N2.r = r;
      particle Nt2(PDG::pdg_proton, PDG::mass_proton); Nt2.r = r;
      e->in.push_back(Nt2); e->out.push_back(N2);
    }
    kaskada k(p, *e, &input);
    k.kaskadaevent(false);
    int npr = 0, nnu = 0, npi = 0; double eprot = 0;
    for (auto &q : e->post) { if (q.pdg == PDG::pdg_proton) { npr++; eprot += q.Ek(); } if (q.pdg == PDG::pdg_neutron) nnu++; if (q.pion()) npi++; }
    printf("%d %d %d %d %d %d %d %d %d %d %.3f\n", e->flag.isTransparent, e->flag.isCorrelated, e->number_of_interactions(),
           npr, nnu, npi, e->nod[0], e->nod[2], e->nod[3], e->nod[8], eprot);
    delete e;
  }
  if (p.kaskada_sudakov) kaskada::print_sud_stats(cerr);
  return 0;
}
