#!/bin/bash
mu=7.93e-07
output_every_Xgen=303
numgen_inN=0.20402
start_output=0.00202
cost=0.0165357839854309
for seed in 100
do
echo "
$seed
$mu
$cost
$output_every_Xgen
$numgen_inN
$start_output
" | ./Code_and_shellscript/HIVevolution_HIV1site >./SimData/Data_T_0.2379_cost_0.0165357839854309.txt
done
