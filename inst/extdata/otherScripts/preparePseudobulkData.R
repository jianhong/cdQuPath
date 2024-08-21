# load the scRNA-seq data
scRNAseqSeu <- readRDS('int.rds')
# do cluster, set higher resolution to get more class
scRNAseqSeu <- FindNeighbors(scRNAseqSeu, dims = 1:20)
scRNAseqSeu <- FindClusters(scRNAseqSeu, resolution = 100)
# create a pseudobulking
pseudo <- AggregateExpression(scRNAseqSeu, assays = "RNA",
                              return.seurat = TRUE,
                              group.by = c("seurat_clusters"))
pseudo <- FindVariableFeatures(pseudo, selection.method = "vst", nfeatures = 200)
genes <- c(VariableFeatures(pseudo), unlist(markers))
pseudo <- subset(pseudo, features=genes)
pseudo$orig.ident <- NULL
saveRDS(pseudo, 'inst/extdata/pseudobulk.scRNAseq.rds')
