#!/bin/bash

directions=("UP" "DOWN")
basedir="/Users/anna_froehlich/ownCloud/single_cell_project/dreamlet/results/ORA/"

for direction in "${directions[@]}"; do
  ${basedir}GO-Figure/gofigure -i ${basedir}input_GO_Figure/mashR_${direction}.tsv -o ${basedir}output_GO_Figure/${direction} --palette cividis --file_type pdf --title mashR_${direction} --outfile_appendix mashR${direction} --description_limit 50
done
