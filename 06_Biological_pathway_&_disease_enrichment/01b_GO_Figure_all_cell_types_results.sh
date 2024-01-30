#!/bin/bash

# upregulated pathways
celltypes_UP=("Astro_PP" "Exc_L2-3" "Exc_L4-6_1" "Exc_L4-6_2" "Exc_L4-6_3" "In_PVALB_Ba" "In_SST" "In_VIP" "Microglia" "Oligodendrocyte" "OPC")
)
basedir="/Users/anna_froehlich/ownCloud/single_cell_project/dreamlet/results/ORA/"

for ct in "${celltypes_UP[@]}"; do
  ${basedir}GO-Figure/gofigure -i ${basedir}input_GO_Figure/${ct}_UP.tsv -o ${basedir}output_GO_Figure/UP --palette cividis --file_type pdf --title UP_${ct} --outfile_appendix ORA_UP_${ct} --description_limit 80
done


# downregulated pathways
celltypes_DOWN=("Astro_FB" "Astro_PP" "Endothelial" "Exc_L2-3" "Exc_L3-5" "Exc_L4-6_1" "Exc_L4-6_2" "Exc_L4-6_3" "Exc_L5-6_1" "In_LAMP5_1" "In_LAMP5_2" "In_PVALB_Ba" "In_PVALB_Ch" "In_RELN" "In_SST" "In_VIP" "Microglia" "Oligodendrocyte" "OPC")

for ct in "${celltypes_DOWN[@]}"; do
  ${basedir}GO-Figure/gofigure -i ${basedir}input_GO_Figure/${ct}_DOWN.tsv -o ${basedir}output_GO_Figure/DOWN --palette cividis --file_type pdf --title DOWN_${ct} --outfile_appendix ORA_DOWN_${ct} --description_limit 70
done