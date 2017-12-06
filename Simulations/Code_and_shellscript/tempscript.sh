#!/bin/bash
mu=1e-05
output_every_Xgen=190
numgen_inN=6.156
start_output=0.076
cost=0.02639
for seed in 100
do
echo "
$seed
$mu
$cost
$output_every_Xgen
$numgen_inN
$start_output
" | ./Code_and_shellscript/HIVevolution_HIV1site5000 >./SimData/Data_T_0.05_cost_0.02639.txt
done
