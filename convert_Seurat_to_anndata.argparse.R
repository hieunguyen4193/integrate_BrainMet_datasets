#####----------------------------------------------------------------------#####
##### 01
#####----------------------------------------------------------------------#####
gc()
rm(list = ls())

library(argparse)
parser <- ArgumentParser()
parser$add_argument("-i", "--input", type="character",
                    help="Path to input seurat object")
parser$add_argument("-o", "--output", type="character",
                    help="Path to save output anndata objects")
parser$add_argument("-r", "--reduction_name", type="character",
                    help="Name of the reduction method")
parser$add_argument("-c", "--cluster_name", type="character",
                    help="Name of the cluster")

args <- parser$parse_args()

path.to.s.obj <- args$input
outputdir <- args$output
reduction.name <- args$reduction_name
cluster.name <- args$cluster_name

s.obj <- readRDS(path.to.s.obj)
PROJECT <- str_replace(basename(path.to.s.obj), ".rds", "")

s.obj$barcode <- colnames(s.obj)

s.obj$UMAP_1 <- s.obj@reductions[[reduction.name]]@cell.embeddings[,1]
s.obj$UMAP_2 <- s.obj@reductions[[reduction.name]]@cell.embeddings[,2]

write.csv(s.obj@reductions[[reduction.name]]@cell.embeddings, 
          file=file.path(outputdir, 
                         sprintf('pca_%s.csv', sprintf("%s_%s_%s", PROJECT, dataset_name, cluster.name))), 
          quote=F, 
          row.names=F)

write.csv(s.obj@meta.data, file=file.path(outputdir, sprintf('metadata_%s.csv', sprintf("%s_%s_%s", PROJECT, dataset_name, cluster.name))), quote=F, row.names=F)

# write expression counts matrix
library(Matrix)
counts_matrix <- GetAssayData(s.obj, assay='SCT', slot='data')
writeMM(counts_matrix, file=file.path(outputdir, sprintf('counts_%s.mtx', sprintf("%s_%s_%s", PROJECT, dataset_name, cluster.name))))

# write gene names
write.table( data.frame('gene'=rownames(counts_matrix)),file=file.path(outputdir, 
                                                                       sprintf('gene_names_%s.csv', 
                                                                               sprintf("%s_%s_%s", PROJECT, dataset_name, cluster.name))), quote=F,row.names=F,col.names=F)
