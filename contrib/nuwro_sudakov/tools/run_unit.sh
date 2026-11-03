#!/bin/bash
# $1 tag, $2 sudakov(0/1), $3 seg, $4 random, $5 headroom, $6 step, $7 seed, $8 sigma
mkdir -p $1 && cd $1
SUD_TEST_SIGMA=$8 SUD_TEST_RANDOM=$4 SUD_TEST_HEADROOM=$5 /home/claude/rr.sh /home/claude/src/nuwro_b/bin/kaskada_test -i /home/claude/src/nuwro_a/data/kaskada.txt -o t.root \
 -p number_of_events=200000 -p beam_particle=2212 -p beam_energy=1438.272 -p random_seed=$7 -p beam_placement=2 \
 -p nucleus_model=0 -p kaskada_NN_corr=0 -p kaskada_sudakov=$2 -p kaskada_sudakov_segment=$3 -p step=$6 > log 2>&1
rm -f t.root
