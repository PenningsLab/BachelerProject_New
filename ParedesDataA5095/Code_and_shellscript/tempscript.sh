#!/bin/bash
mu=7e-07
output_every_Xgen=343
numgen_inN=0.230953333333333
start_output=0.00228666666666667
cost=0.0145965306302669
for seed in 100
do
echo "
$seed
$mu
$cost
$output_every_Xgen
$numgen_inN
$start_output
" | ./Code_and_shellscript/HIVevolution_HIV1site >./SimData/Data_T_0.3172776_cost_0.0145965306302669.txt
done
