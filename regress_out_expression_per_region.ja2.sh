#!/bin/bash
#$ -V
#$ -cwd
#$ -pe shared 1
#$ -l h_data=80G,h_rt=6:00:00
#$ -m a
#$ -o regress_out_expression_per_region.out
#$ -e regress_out_expression_per_region.err

ja=regress_out_expression_per_region.job
log_file_prefix="regress_out_expression_per_region."

module load /u/local/Modules/modulefiles/R/4.2.2

PARMS=($(awk "NR==$SGE_TASK_ID" $ja))
infile=${PARMS[0]}
outfile=${PARMS[1]}

Rscript regress_out_expression_per_region.ROSMAP.R ${infile} ${outfile} 1>log/${log_file_prefix}${SGE_TASK_ID}.out 2>log/${log_file_prefix}${SGE_TASK_ID}.err

echo "job completed :D" >> log/${log_file_prefix}${SGE_TASK_ID}.out

#regress_out_expression_per_region.ja2.sh

