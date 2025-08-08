#creates job array

ja=regress_out_expression_per_region.job
inpath=/u/project/gxxiao/tcspence/projects/ediTWAS/data/ROSMAP/bulk_RNA_seq/phenotype_tables/region_based/rosmap_elaine_editing_sites.no_race_in_model.regression_all_together._
file_suffix="_phenotype_pcs_regressed_out.txt" #everything up till ".txt"

fixed_inpath=/u/project/gxxiao/tcspence/projects/ediTWAS/data/ROSMAP/bulk_RNA_seq/phenotype_tables/region_based/rosmap_elaine_editing_sites.regression_all_together_

################DONE CHANGING PARAMETERS##################
rm $ja #remove old job file

pcs=(1 2 3 4 5)
#pcs=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20)

for num_pcs in ${pcs[*]}
do
	infile=$(ls ${inpath}${num_pcs}${file_suffix})
	echo "${infile} ${fixed_inpath}${num_pcs}${file_suffix}.expression_regressed_per_region.txt" >> $ja
done

echo job completed

