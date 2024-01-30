#!/bin/bash

celltypes=("Astro_FB" "Astro_PP" "Endothelial" "Exc_L2-3" "Exc_L3-5" "Exc_L4-6_1" "Exc_L4-6_2" "Exc_L4-6_3" "Exc_L5-6_1" "Exc_L5-6_2" "Exc_L5-6_HTR2C" "In_LAMP5_1" "In_LAMP5_2" "In_PVALB_Ba" "In_PVALB_Ch" "In_RELN" "In_SST" "In_VIP" "Microglia" "Oligodendrocyte" "OPC")
#traits=("hypertension")
#ct=("Oligodendrocyte")
traits=("BIP2021" "MDD" "SCZ" "AD" "hypertension")
basedir="/Users/anna_froehlich/ownCloud/single_cell_project/H-MAGMA/"

for trait in "${traits[@]}"; do
  for ct in "${celltypes[@]}"; do
    ${basedir}magma_v1/magma --gene-results ${basedir}SNP_to_gene_annotation/adult_brain_${trait}.genes.raw --gene-covar ${basedir}MAGMAinput/Age/${ct}_gene_covar.txt --out ${basedir}MAGMAoutput/Age/${trait}_${ct}_magmaEnrichment
  done
done
