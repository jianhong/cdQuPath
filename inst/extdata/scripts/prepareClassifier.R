d <- read.csv('~/Downloads/Wu_EMBO_countr_matrix.csv')
meta <- read.csv('~/Downloads/Wu_EMBO_metadata.csv', skip=2, header=FALSE)
celltypes <- meta$V7
identical(colnames(d)[-1], meta$V1) # TRUE
genes <- readLines('~/Downloads/Wu_EMBOJ_genes.tsv')
identical(genes, meta$V1) # FALSE
markers <- readRDS(system.files('extdata', 'markers.name.map.rds'), package='cdQuPath')
m <- unlist(markers)
rownames(d) <- d[, 1]
m <- m[m %in% rownames(d)]
dd <- d[m, ]
dd <- dd[, -1]
m1 <- lapply(markers[sapply(markers, function(.ele) any(.ele %in% rownames(dd)))], function(.ele){
  colSums(dd[.ele[.ele %in% rownames(dd)], , drop=FALSE])
})
m1[1:5]
m1 <- do.call(rbind, m1)
rowSums(m1)
seu <- CreateSeuratObject(m1)
seu <- fitGMM(seu)
gmm <- GetAssayData(seu, assay = 'GMM', layer = 'data')
identical(colnames(gmm), meta$V1)
classifier <- buildClassifier(gmm, features = rownames(gmm), celltypes = celltypes)
saveRDS(classifier, 'inst/extdata/classifier.SCP1106.rds')

seu <- readRDS('/Users/ouj/Library/CloudStorage/Box-Box/Reanalysis_filteredcellsremoved/resolution 0.3/int.rds')
d <- GetAssayData(seu, assay = 'RNA', layer = 'count')
dim(d)
dd <- d[m, ]
dim(dd)
m1 <- lapply(markers[sapply(markers, function(.ele) any(.ele %in% rownames(dd)))], function(.ele){
  colSums(dd[.ele[.ele %in% rownames(dd)], , drop=FALSE])
})
m1 <- do.call(rbind, m1)
m1 <- log2(m1+1)
rowSums(m1)
seu1 <- CreateSeuratObject(m1)
seu1 <- fitGMM(seu1)
gmm <- GetAssayData(seu1, assay = 'GMM', layer = 'data')
identical(colnames(gmm), colnames(seu))
clusters <- seu$seurat_clusters
celltypes <- c('0'='Fibroblastic',
               '1'="Chondroblastic",
               '2'='Myeloid',
               '3'='Chondroblastic',
               '4'='Chondroblastic',
               '5'='Pericytes',
               '6'='Osteoblastic',
               '7'='Proliferating',
               '8'='Endothelial',
               '9'='Osteoclastic',
               '10'='T/NK cells',
               '11'='Mast cells')
cts <- celltypes[as.character(clusters)]
gmm[is.na(gmm)] <- 0
classifier <- buildClassifier(gmm, features = rownames(gmm), celltypes = cts)
saveRDS(classifier, 'inst/extdata/classifier.CSA.res0.3.rds')
