# regress-expression-per-region

# AEI Expression Regression Script

This R script performs **parallelized regression of gene expression levels from RNA editing data**, specifically region-based Alu Editing Index (AEI) values. The script was used on ROSMAP bulk RNA-seq datasets aligned to the GRCh38 genome and is designed to prepare RNA editing phenotypes for downstream association testing (e.g., edQTL or ediTWAS analyses).

## Background

RNA editing levels can be confounded by gene expression, as highly expressed genes tend to produce more reads and therefore more detectable editing. This can introduce spurious correlations between genotype and editing that are actually driven by differences in gene expression. To control for this, the script regresses out gene expression from each Alu element’s AEI value. Residuals from this regression represent editing levels independent of expression, and are used as the phenotype in downstream genetic association analyses.

## Inputs

- **editing_file**: A TSV with AEI values (Alu Editing Index) per region (Alu element) per sample.
- **expression_matrix_file**: A large, TPM-normalized gene expression matrix from ROSMAP bulk RNA-seq data.
- **ensembl_annotated_alu_file**: A BED file mapping Alu regions to host genes via Ensembl IDs (created using BEDtools intersect).
- **covariate and conversion tables**: Used to align metadata and match sample names across input datasets.

## Functions

- Cleans and processes the expression matrix (e.g., NA filtering, imputation, Ensembl ID matching).
- Merges editing data with Alu–host gene annotations.
- Matches gene expression rows by Ensembl ID and ensures sample order consistency.
- Runs linear regression for each Alu element:
  ```
  AEI_region ~ Expression_gene
  ```
  Extracts residuals as the deconfounded editing signal.
- Uses `foreach` and `doParallel` for parallel processing across all Alu elements.

**The ja1/ja2 job array system (as demonstrated in the ".sh" files) parallelizes these functions by sample, allowing you to run this for multiple Alus in multiple samples all at once, greatly reducing runtime.**

## Output

A tab-separated file where:
- Rows = Alu regions
- Columns = Samples
- Values = Expression-residualized AEI per sample per Alu region

This output was used as the **phenotype input** for tensorQTL mapping.
