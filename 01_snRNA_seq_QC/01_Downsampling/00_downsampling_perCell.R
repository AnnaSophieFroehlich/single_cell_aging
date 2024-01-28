##################################################
## Project: Single nuclei postmortem brain
## Date: 02.08.2021
## Author: Nathalie
##################################################
# Downsampling of reads per cell to 75% quantile
# of mean reads per cell 

# run with conda environment downsampling

library(DropletUtils)
library(data.table)
library(dplyr)

basepath_batch1 <- "/psycl/g/mpsngs/HiSeq_Helmholtz/20201029_Anna_Froehlich_10X_RNAseq/01_cellranger_v6/"
basepath_batch2 <- "/psycl/g/mpsngs/HiSeq_Helmholtz/20210324_Anna_Froehlich_10X_RNAseq/01_cellranger_v6/"

# read file with read proportions
prop <-
  read.table(
    "/psycl/g/mpsngs/HiSeq_Helmholtz/20210324_Anna_Froehlich_10X_RNAseq/01_cellranger_v6/downsample_prop.csv",
    header = TRUE
  )

# target number of reads --> 75% quantile
target.reads <- 14786.00


# downsample the reads from molecule file
for (i in 1:nrow(prop)){
  
  output_dir <- file.path(
    "/psycl/g/mpsngs/HiSeq_Helmholtz/20210324_Anna_Froehlich_10X_RNAseq/03_downsampled/",
    prop$Sample[i]
  )
  
  # create folder for output if doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  # path to input file on which downsampling happens and feature file
  if (prop$Main.Batch[i] == 1){
    mol.info.file <- file.path(basepath_batch1, prop$Sample[i], "outs/molecule_info.h5")
    feature.file <- file.path(basepath_batch1, prop$Sample[i], "outs/raw_feature_bc_matrix/features.tsv.gz")
  } else {
    mol.info.file <- file.path(basepath_batch2, prop$Sample[i], "outs/molecule_info.h5")
    feature.file <- file.path(basepath_batch2, prop$Sample[i], "outs/raw_feature_bc_matrix/features.tsv.gz")
  }
  
  # read mol info file
  mol.info <- read10xMolInfo(mol.info.file)
  str(mol.info)
  sum(mol.info$data$reads)
  
  # get reads per cell
  reads.per.cell <- as.data.frame(mol.info$data) %>%
    group_by(cell) %>%
    summarise(reads.cell = sum(reads))
  
  # add column with downsampling proportion per cell
  reads.per.cell <- reads.per.cell %>%
    mutate(downsample.prop = target.reads/reads.per.cell$reads.cell)
  reads.per.cell$downsample.prop[reads.per.cell$downsample.prop >= 1] <- 1
  
  # read feature file
  features <- fread(feature.file, header = FALSE)
  
  # UMI counts per cell
  umi_cell <- sort(table(mol.info$data$cell), decreasing = TRUE)
  head(umi_cell)
  q99 <- quantile(umi_cell[1:8000], 0.99)
  cutoff <- q99/10
  cutoff
  
  # downsample reads
  data.sampling <- downsampleReads(mol.info.file, 
                                   prop = reads.per.cell$downsample.prop,
                                   bycol = TRUE)
  
  # barcode ranks
  br.out <- barcodeRanks(data.sampling)
  str(br.out)
  
  # Making a plot.
  png(file.path(output_dir, "rank_barcode_perCell.png"), width = 500, height = 500)
  plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
  o <- order(br.out$rank)
  lines(br.out$rank[o], br.out$fitted[o], col="red")
  
  abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
  abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
  legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
         legend=c("knee", "inflection"))
  dev.off()
  
  # don't take mitochondrial and ribosomal genes into account for emptyDrops
  keep <- !startsWith(features$V2, "MT-") & !startsWith(features$V2, "RP")
  
  # Identify likely cell-containing droplets.
  #out <- emptyDrops(data.sampling)    
  #out <- emptyDrops(data.sampling[keep,])    
  #out <- emptyDrops(data.sampling, retain = Inf)    
  #out <- emptyDrops(data.sampling, retain = cutoff)   
  out <- emptyDrops(data.sampling[keep,], retain = cutoff)    
  #out <- emptyDrops(data.sampling, retain = cutoff, lower = 50) 
  out
  is.cell <- out$FDR <= 0.01
  sum(is.cell, na.rm=TRUE)
  # Subsetting the matrix to the cell-containing droplets.
  # (using 'which()' to handle NAs smoothly).
  cell.counts <- data.sampling[,which(is.cell),drop=FALSE]
  dim(cell.counts)
  # Check if p-values are lower-bounded by 'niters'
  # (increase 'niters' if any Limited==TRUE and Sig==FALSE)
  table(Sig=is.cell, Limited=out$Limited)
  
  # diagnostic plot for emptyDrops
  png(file.path(output_dir, "emptyDrops_diagnostics_perCell.png"), width = 700, height = 500)
  plot(out$Total, -out$LogProb, col=ifelse(is.cell, "red", "black"),
       xlab="Total UMI count", ylab="-Log Probability")
  dev.off()
  
  # Write 10X counts to file
  write10xCounts(
    file.path(output_dir, "bc_feature_matrix_downsampling_perCell.h5"),
    cell.counts,
    gene.symbol = features$V2,
    overwrite = TRUE,
    type = "HDF5"
  )
  
}

