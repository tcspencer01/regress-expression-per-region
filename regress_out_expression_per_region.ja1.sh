#creates job array

ja=regress_out_expression_per_region.job
inpath=editing_sites.regression_all_together_
file_suffix="_phenotype_pcs_regressed_out.txt" #everything up till ".txt"

################DONE CHANGING PARAMETERS##################
rm $ja #remove old job file

pcs=(1 2 3 4 5)

for num_pcs in ${pcs[*]}
do
	infile=$(ls ${inpath}${num_pcs}${file_suffix})
	echo "${infile} ${inpath}${num_pcs}${file_suffix}.expression_regressed_per_region.txt" >> $ja
done

echo job completed

