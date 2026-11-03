#!/usr/bin/env bash
# Runs the Sudakov validation suite (Step 4 a-e). Usage: run_all.sh <build dir> <out dir> <qmc_configs.txt> <qmc_density.txt>
set -euo pipefail
B="$(cd "$1" && pwd)"; O="$(mkdir -p "$2" && cd "$2" && pwd)"; CFG="$(realpath "$3")"; DEN="$(realpath "$4")"
cd "$B"; V="$B/src/Achilles/achilles-sudakov-validate"
mkdir -p "$O"; N=200000
run() { local name=$1; shift; "$V" "$@" out="$O/$name.csv" 2> >(grep "events=" > "$O/$name.diag") >/dev/null; }
# (a) closed forms, Track A (headroom exercises the veto)
run a_ext_A sphere-ext track=A sigma=40 headroom=1.5 nev=$N seed=11
run a_int_A sphere-int track=A sigma=40 headroom=1.5 nev=$N seed=12
run a_ext_A20 sphere-ext track=A sigma=20 headroom=1.0 nev=$N seed=13
run a_ext_A30 sphere-ext track=A sigma=30 headroom=1.0 nev=$N seed=14
# (a) Track B on Poisson configurations: N|b ~ Poisson((T*p)(b)) exactly, any profile
for prof in Atom Box Gauss; do
  run a_ext_B_k1_$prof sphere-ext track=B sigma=40 law=Kappa1 profile=$prof nev=$N seed=21
  run a_ext_B_eik_$prof sphere-ext track=B sigma=40 law=Eikonal x0=2.6 profile=$prof headroom=1.3 nev=$N seed=22
done
# hazard thinning (Pauli-like constant factor a=0.5)
run a_int_A_pauli sphere-int track=A sigma=40 pauli=0.5 nev=$N seed=31
run a_ext_B_pauli sphere-ext track=B sigma=40 law=Eikonal x0=2.6 profile=Box pauli=0.5 nev=$N seed=32
# (b) toy enumeration
for t in 2 3 4; do for prof in Atom Box Gauss; do
  run b_toy${t}_$prof toy toy=$t sigma=40 law=Kappa1 profile=$prof nev=$N seed=4$t
done; done
for prof in Atom Box Gauss; do run b_toy3int_$prof toy toy=3 internal=1 sigma=40 law=Kappa1 profile=$prof nev=$N seed=49; done
# (c) violated majorant
run c_ext_A_h08 sphere-ext track=A sigma=40 headroom=0.8 nev=$N seed=51
run c_int_A_h08 sphere-int track=A sigma=40 headroom=0.8 nev=$N seed=52
run c_ext_A_rho016 qmc-ext track=A sigma=40 rhohat=0.16 density="$DEN" configs="$CFG" nev=$N seed=53
run c_int_A_rho016 qmc-int track=A sigma=40 rhohat=0.16 density="$DEN" configs="$CFG" nev=$N seed=54
run c_ext_B_h07 sphere-ext track=B sigma=40 law=Eikonal profile=Box headroom=0.7 nev=$N seed=55
# (d) reweighting closure
for c in 1.25 2 4; do run d_int_A_c$c sphere-int track=A sigma=40 headroom=$c vars=0.8,1.2 nev=$N seed=61; done
run d_ext_B_rw sphere-ext track=B sigma=40 law=Eikonal profile=Gauss headroom=1.0 vars=0.8,1.2 nev=$N seed=62
run d_ext_B_32 sphere-ext track=B sigma=32 law=Eikonal profile=Gauss nev=$N seed=63
run d_ext_B_48 sphere-ext track=B sigma=48 law=Eikonal profile=Gauss nev=$N seed=64
run d_int_A_pvar sphere-int track=A sigma=40 pauli=0.6 pvars=0.3,0.9 nev=$N seed=65
# (e) A vs B on QMC (sigma = 40 mb)
W=$(python3 -c "import math;print(math.sqrt(4/(2*math.pi)))")
run e_ext_A_rho qmc-ext track=A sigma=40 density="$DEN" configs="$CFG" nev=$N seed=71
run e_ext_A_fold qmc-ext track=A sigma=40 fold=$W density="$DEN" configs="$CFG" nev=$N seed=72
run e_ext_A_foldFA qmc-ext track=A sigma=40 fold=$W finiteA=1 density="$DEN" configs="$CFG" nev=$N seed=73
for dec in 0 1; do
  run e_ext_B_k1_dec$dec qmc-ext track=B sigma=40 law=Kappa1 profile=Box decorrelate=$dec configs="$CFG" nev=$N seed=74
  run e_ext_B_eik_dec$dec qmc-ext track=B sigma=40 law=Eikonal x0=2.6 profile=Box headroom=1.3 decorrelate=$dec configs="$CFG" nev=$N seed=75
  for prof in Atom Box Gauss; do
    run e_int_B_${prof}_dec$dec qmc-int track=B sigma=40 law=Eikonal x0=2.6 profile=$prof headroom=1.3 decorrelate=$dec configs="$CFG" nev=$N seed=76
  done
done
run e_int_A_rho qmc-int track=A sigma=40 density="$DEN" configs="$CFG" nev=$N seed=77
run e_int_A_foldFA qmc-int track=A sigma=40 fold=$W finiteA=1 density="$DEN" configs="$CFG" nev=$N seed=78
echo done
