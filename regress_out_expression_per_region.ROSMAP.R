#Regress expression out of phenotype (editing) matrix
#ROSMAP hg38 (region-based)

rm(list=ls())

#VERY IMPORTANT: This needs to be run on hoffman2 because the expression file is so huge

library(missMethods)
library(tidyr)
library(dplyr)

args = commandArgs(trailingOnly = TRUE)
editing_file = args[1]
outfile = args[2]

# ####################TO TEST ON HOFFMAN2 INTERACTIVE NODE (48G)
# editing_file = "/u/project/gxxiao/gxxiao2/tcspence/ROSMAP/grch38/region_based/output/rosmap_elaine_editing_sites.no_race_in_model.regression_all_together._1_phenotype_pcs_regressed_out.txt"
# outfile = "/u/project/gxxiao/gxxiao2/tcspence/ROSMAP/grch38/region_based/output/TEST.txt"
# ##############################################################

ensembl_annotated_alu_file = "/u/project/gxxiao/tcspence/projects/ediTWAS/annotation/hg38.fa.Alus.RNAEditingIndexer.stranded.ensembl_host_gene_overlapped_less_stringent.bed"

expression_matrix_file = "/u/project/gxxiao/tcspence/projects/ediTWAS/data/ROSMAP/bulk_RNA_seq/phenotype_tables/expression/TPM_normalized_expression.ROSMAP_bulk.552_samples.expanded_info.tsv"

bam_conversion_file = "/u/project/gxxiao/tcspence/projects/ediTWAS/dataset_info/ROSMAP/bam_to_rosmap_id_conversion_table.txt"

conversion_file = "/u/project/gxxiao/tcspence/projects/ediTWAS/dataset_info/ROSMAP/ROSMAP_biospecimen_metadata.ID_conversion.txt"

vcf_header_file = "/u/project/gxxiao/tcspence/projects/ediTWAS/dataset_info/ROSMAP/rosmap_vcf_header.txt"

full_covariates_file = "/u/project/gxxiao/tcspence/projects/ediTWAS/dataset_info/ROSMAP/known_covariates_site_based_tensorQTL_input.txt"

######################################################################

input_editing = read.table(editing_file,sep="\t",header=TRUE,stringsAsFactors=FALSE,quote="\"",check.names=FALSE)
ensembl_annotated_alu_table = read.table(ensembl_annotated_alu_file,sep="\t",header=FALSE,stringsAsFactors=FALSE,quote="\"",check.names=FALSE)
colnames(ensembl_annotated_alu_table) <- c("alu_chr", "alu_start", "alu_end", "alu_score", "alu_name", "alu_strand", "gene_chr", "gene_start", "gene_end", "gene_id", "gene_family", "gene_strand")

num_samples = ncol(input_editing)-1

#Load known covariates
full_covariates = read.table(full_covariates_file,sep="\t",header=TRUE,check.names=FALSE)
conversion_table = read.table(conversion_file,sep="\t",header=TRUE)
bam_conversion_table = read.table(bam_conversion_file,sep="\t",header=TRUE,stringsAsFactors=FALSE,quote="\"",check.names=FALSE)

bam_conversion_table = bam_conversion_table[!grepl("DLPFC", bam_conversion_table$name),]
bam_conversion_table$name = gsub('.{4}$', '', bam_conversion_table$name)

conversion_table = subset(conversion_table, specimenID %in% colnames(full_covariates))
conversion_table = merge(conversion_table, bam_conversion_table, by="individualID")
colnames(conversion_table)[3] = "bam_name"

expression_matrix = read.table(expression_matrix_file,sep="\t",header=TRUE,check.names=FALSE)
expression_matrix_backup = expression_matrix
#expression_matrix = expression_matrix_backup
print(paste("Raw expression matrix nrow:", nrow(expression_matrix)))

expression_matrix = subset(expression_matrix, !(is.na(Ensembl_Gene_ID)))
print(paste("Raw expression matrix nrow after filtering NA ensembl IDs:", nrow(expression_matrix)))

expression_matrix = expression_matrix[rowSums(is.na(expression_matrix[,1:(ncol(expression_matrix)-8)])) != (ncol(expression_matrix)-8), ]
print(paste("NA-filtered expression matrix nrow:", nrow(expression_matrix)))

expression_matrix_info_cols = expression_matrix[,(ncol(expression_matrix)-7):(ncol(expression_matrix))]

expression_matrix = expression_matrix[rowSums(is.na(expression_matrix)) != (ncol(expression_matrix)), ]
print(paste("NA-filtered expression matrix nrow (after filtering for sample names):", nrow(expression_matrix)))

expression_matrix = expression_matrix[,order(colnames(expression_matrix))]

#Make the matrix numeric
cols = colnames(expression_matrix)
rows = rownames(expression_matrix)
expression_matrix = data.frame(lapply(expression_matrix, function(x) as.numeric(as.character(x))))
colnames(expression_matrix) = cols
rownames(expression_matrix) = rows

expression_matrix = impute_mean(expression_matrix, type="rowwise") #Impute expression for samples who don't have expression measured for a particular gene

#any(is.na(expression_matrix)) #debugging

expression_matrix = cbind(expression_matrix_info_cols$Ensembl_Gene_ID, expression_matrix)
colnames(expression_matrix)[1] = "ensembl_gene_id"

if(length(unique(expression_matrix$ensembl_gene_id)) == nrow(expression_matrix)){
  print("No ensembl ids appear twice in the input expression matrix - yay!")
} else{
  print("Duplicate ensembl ids in input expression matrix --> Could cause problems")
  stop("ERROR: Duplicate ensembl ids")
}

#which(duplicated((expression_matrix$ensembl_gene_id))) #debugging

gc() #Since the merge coming up takes a lot of memory

input_editing$gene_id = separate_wider_delim(input_editing, region, ":", names=c("chr", "pos", "alu_strand", "alu_type", "gene_id", "gene_strand"))$gene_id
ensembl_thing = distinct(ensembl_annotated_alu_table[,7:12])
input_editing = merge(input_editing, ensembl_thing, by="gene_id")
#input_editing = distinct(input_editing) #in case it duplicated
input_editing = subset(input_editing, gene_id %in% expression_matrix$ensembl_gene_id)

input_editing_regionInfo = input_editing[,1:2]
input_editing = input_editing[3:ncol(input_editing)]
input_editing = input_editing[,order(colnames(input_editing))]
input_editing = cbind(input_editing_regionInfo, input_editing)

#expression_matrix_backup2 = expression_matrix
expression_matrix = expression_matrix[, colnames(expression_matrix) %in% c("ensembl_gene_id", colnames(input_editing))]

# #DEBUGGING
# length(colnames(input_editing)[2:(ncol(input_editing)-6)])
# col1 = colnames(input_editing)[2:(ncol(input_editing)-6)]
# col2 = colnames(expression_matrix)[2:ncol(expression_matrix)]
# for(i in c(1:length(col1))){ #Shows the colnames are the same
#   if(col1[i] != col2[i]){
#     print(paste(i, col1[i], col2[i]))
#   }
# }

#PARALLELIZED VERSION:
#Based on: https://stackoverflow.com/questions/38318139/run-a-for-loop-in-parallel-in-r
library(foreach)
library(doParallel)

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

#save some variables to reduce runtime
expression_matrix_ncol = ncol(expression_matrix)
input_editing_ncol = ncol(input_editing)

#now run parallel for loop
finalMatrix <- foreach(i=1:nrow(input_editing), .combine=rbind) %dopar% {
#finalMatrix <- foreach(i=1:3, .combine=rbind) %dopar% { #debugging
  
  current_ensembl_id = input_editing$gene_id[i]
  
  linear_model = lm(as.numeric(as.vector(t(input_editing[i,8:(input_editing_ncol)]))) ~ as.numeric(as.vector(t(subset(expression_matrix, ensembl_gene_id==current_ensembl_id)[,2:expression_matrix_ncol]))), na.action="na.exclude") #Regress expression vector out of editing vector for this region --> This assumes each ensembl id only appears once
  
  linear_model_residuals = linear_model$residuals
  linear_model_residuals #equivalent to rbind with final matrix
  
}
#The above parallelized version preserves row order, so you can just move the old col names and then cbind regions back in:
finalMatrix = cbind(input_editing[,2], finalMatrix)
colnames(finalMatrix) = c("region", colnames(input_editing)[8:(length(input_editing))])

write.table(finalMatrix, file = outfile, sep = "\t", quote=FALSE, row.names = F, col.names = T)

print("Job completed in R :D")

