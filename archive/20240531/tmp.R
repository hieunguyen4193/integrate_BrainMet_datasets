##### tmp
s.obj.merge <- readRDS("/media/hieunguyen/HNHD01/data/UKK/UKK_draft/Lung/integrated_dataset/merge_myloid_from_integrated_dataset_and_control_dataset.rds")
s.obj.merge <- JoinLayers(s.obj.merge)
tmp.s.obj.myeloid <- readRDS("/media/hieunguyen/HNHD01/data/UKK/UKK_draft/Lung/integrated_dataset/integrated_datasets_Myeloid_only.rds")

p.myeloid <- DimPlot(object = tmp.s.obj.myeloid, reduction = "umap", label = TRUE, label.box = TRUE, repel = TRUE)
p.final <- DimPlot(object = s.obj.merge, reduction = "RNA_UMAP", group.by = "final.celltype", label = TRUE, label.box = TRUE, repel = TRUE)

ggsave(plot = p.myeloid, filename = "myeloid_cells_from_integrated_dataset.svg", path = path.to.03.output, device = "svg", width = 14, height = 10, dpi = 300)
ggsave(plot = p.final, filename = "final.svg", path = path.to.03.output, device = "svg", width = 14, height = 10, dpi = 300)

kinase.genelist <- read.csv("/media/hieunguyen/HNHD01/data/UKK/geneSet/TK_genes.csv")$TK_gene
kinase.genelist <- intersect(kinase.genelist, row.names(s.obj.merge))
check.celltypes <- c("Proliferating Mac", 
                     "Alveolar Mac", 
                     "cDC2/moDCs",
                     "Monocytes",
                     "Low quality Mac", 
                     "Neutrophils",
                     "Lipid-associated Mac")

all.s.obj <- hash()
for (input.celltype in check.celltypes){
  all.s.obj[[input.celltype]] <- subset(s.obj.merge, final.celltype == input.celltype)
}

for (input.celltype in check.celltypes){
  print(sprintf("working on %s", input.celltype))
  dir.create(file.path(path.to.03.output, "violin_plot", input.celltype), showWarnings = FALSE, recursive = TRUE)
  for (gene.id in kinase.genelist){
    tmp <-  all.s.obj[[input.celltype]]
    counts <- GetAssayData(tmp, assay = "RNA", slot = "counts")[gene.id, ]
    print(sprintf("working on gene %s", gene.id))
    keep.cells <- names(counts[counts > 0])
    if (length(keep.cells) != 0){
      if (file.exists(file.path(path.to.03.output, "violin_plot", input.celltype, sprintf("violinplot_%s_%s.png", gene.id, str_replace(input.celltype, "/", "-"))))== FALSE){
        p <- VlnPlot(object = subset(tmp, cells = keep.cells), features = c(gene.id), group.by = "sampleid", pt.size = 1) +
          stat_compare_means(method = "wilcox")
        ggsave(plot = p, filename = sprintf("violinplot_%s_%s.png", gene.id, str_replace(input.celltype, "/", "-")), 
               path = file.path(path.to.03.output, "violin_plot", input.celltype), device = "png", width = 14, height = 10, dpi = 300)          
      }
    }
  }
}
new.gene.list <- read.csv(file.path("/media/hieunguyen/HNHD01/data/UKK/geneSet/new_gene_set_20240409.csv"))$Gene

for (input.celltype in check.celltypes){
  print(sprintf("working on %s", input.celltype))
  dir.create(file.path(path.to.03.output, "violin_plot", input.celltype), showWarnings = FALSE, recursive = TRUE)
  for (gene.id in intersect(new.gene.list, row.names(s.obj.merge))){
    tmp <-  all.s.obj[[input.celltype]]
    counts <- GetAssayData(tmp, assay = "RNA", slot = "counts")[gene.id, ]
    print(sprintf("working on gene %s", gene.id))
    keep.cells <- names(counts[counts > 0])
    if (length(keep.cells) > 1){
      if (file.exists(file.path(path.to.03.output, "violin_plot", input.celltype, sprintf("violinplot_%s_%s.png", gene.id, str_replace(input.celltype, "/", "-"))))== FALSE){
        p <- VlnPlot(object = subset(tmp, cells = keep.cells), features = c(gene.id), group.by = "sampleid", pt.size = 1) +
          stat_compare_means(method = "wilcox")
        ggsave(plot = p, filename = sprintf("violinplot_%s_%s.png", gene.id, str_replace(input.celltype, "/", "-")), 
               path = file.path(path.to.03.output, "violin_plot", input.celltype), device = "png", width = 14, height = 10, dpi = 300)          
      }
    }
  }
}