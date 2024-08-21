markers <- readRDS(system.file('extdata', 'markers.name.map.rds',
                                package='cdQuPath'))

seu <- readRDS('/Users/ouj/Library/CloudStorage/Box-Box/Reanalysis_filteredcellsremoved/resolution 0.3/int.rds')
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
seu$celltype <- as.character(cts)
pseudo <- AggregateExpression(seu, assays = 'RNA', return.seurat = TRUE, group.by="celltype")
d <- GetAssayData(pseudo, assay = 'RNA', layer = 'count')
dim(d)
dd <- d[c(m, 'IBSP'), ]
dim(dd)
dd1 <- as.matrix(dd)
dd1 <- t(t(dd1)/colSums(dd1))
colSums(dd1)
library(ggplot2)
library(reshape2)
dd1 <- melt(dd1)
colnames(dd1) <- c('marker', 'celltype', 'value')
ggplot(dd1, aes(fill=marker, y=value, x=celltype)) + 
  geom_bar(position="stack", stat="identity")
ggsave('pattern.withIBSP.pdf')
dd <- d[m, ]
dd1 <- as.matrix(dd)
dd1 <- t(t(dd1)/colSums(dd1))
colSums(dd1)
library(ggplot2)
library(reshape2)
dd1 <- melt(dd1)
colnames(dd1) <- c('marker', 'celltype', 'percentage')
ggplot(dd1, aes(fill=marker, y=percentage, x=celltype)) + 
  geom_bar(position="stack", stat="identity")
ggsave('pattern.withoutIBSP.pdf')

m1 <- lapply(markers[sapply(markers, function(.ele) any(.ele %in% rownames(dd)))], function(.ele){
  colSums(dd[.ele[.ele %in% rownames(dd)], , drop=FALSE])
})
m1 <- do.call(rbind, m1)
rowSums(m1)
seu1 <- CreateSeuratObject(m1)
seu1 <- fitGMM(seu1)
gmm <- GetAssayData(seu1, assay = 'GMM', layer = 'data')
gmm
gmm[is.na(gmm)] <- 0
classifier <- buildClassifier(gmm, features = rownames(gmm), celltypes = colnames(gmm))
classifier$y

saveRDS(classifier, 'inst/extdata/classifier.CSA.res0.3.rds')

