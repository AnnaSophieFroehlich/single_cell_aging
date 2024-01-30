#!/bin/bash

# Gene analysis - SNP p-values
# Map GWAS hits to genes
# get genes.raw file

basedir="/Users/anna_froehlich/ownCloud/single_cell_project/H-MAGMA/"

# Bipolar disorder
#N_GWAS= 413466

${basedir}magma_v1/magma --bfile ${basedir}g1000_eur/g1000_eur --pval ${basedir}GWAS_summary_stats/BIP2021.sumstats N=413466 --gene-annot ${basedir}SNP_to_gene_annotation/Adult_brain.genes.annot --out ${basedir}SNP_to_gene_annotation/adult_brain_BIP2021


# Major depressive disorder
# N_GWAS= 500199

${basedir}magma_v1/magma --bfile ${basedir}g1000_eur/g1000_eur --pval ${basedir}GWAS_summary_stats/MDD.sumstats N=500199 --gene-annot ${basedir}SNP_to_gene_annotation/Adult_brain.genes.annot --out ${basedir}SNP_to_gene_annotation/adult_brain_MDD


# Schizophrenia
# N_GWAS= 127906

${basedir}magma_v1/magma --bfile ${basedir}g1000_eur/g1000_eur --pval ${basedir}GWAS_summary_stats/SCZ2022.sumstats N=127906 --gene-annot ${basedir}SNP_to_gene_annotation/Adult_brain.genes.annot --out ${basedir}SNP_to_gene_annotation/adult_brain_SCZ


# Alzheimer's disease
# N_GWAS= 455258

${basedir}magma_v1/magma --bfile ${basedir}g1000_eur/g1000_eur --pval ${basedir}GWAS_summary_stats/AD.sumstats N=455258 --gene-annot ${basedir}SNP_to_gene_annotation/Adult_brain.genes.annot --out ${basedir}SNP_to_gene_annotation/adult_brain_AD


# Hypertension
# N_GWAS= 361194

${basedir}magma_v1/magma --bfile ${basedir}g1000_eur/g1000_eur --pval ${basedir}GWAS_summary_stats/hypertension.sumstats N=361194 --gene-annot ${basedir}SNP_to_gene_annotation/Adult_brain.genes.annot --out ${basedir}SNP_to_gene_annotation/adult_brain_hypertension



