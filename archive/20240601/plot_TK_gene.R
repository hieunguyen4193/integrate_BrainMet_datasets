library(Seurat)

path.to.s.obj <- "/media/hieunguyen/HNHD01/outdir/BrainMet_rerun_without_SCT_20240615/BrainMet_20240615_SeuratV5_GSE186344_GSM5645888/GSE186344_GSM5645888/s8a_output/BrainMet_20240615_SeuratV5_GSE186344_GSM5645888.output.s8a.rds"
path.to.s.obj.sct <- "/media/hieunguyen/HNHD01/outdir/BrainMet_SeuratV5_GSE186344_GSM5645888/GSE186344_GSM5645888/s8a_output/BrainMet_SeuratV5_GSE186344_GSM5645888.output.s8a.rds"
tk.gene.set <- read.csv("/media/hieunguyen/HNHD01/data/UKK/geneSet/TK_genes.csv")$TK_gene
s.obj <- readRDS(path.to.s.obj)
s.obj.sct <- readRDS(path.to.s.obj.sct)


p <- FeaturePlot(object = s.obj, reduction = "RNA_UMAP", features = head(tk.gene.set, 9))
p.sct <- FeaturePlot(object = s.obj.sct, reduction = "RNA_UMAP", features = head(tk.gene.set, 9))
