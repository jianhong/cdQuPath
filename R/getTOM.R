#' Topological overlap matrix similarity
#' @description
#' Calculate TOM by WGCNA
#' 
#' @param seu A Seurat object
#' @param softPower Soft thresholding power.
#' @param assay The expression assay name. Do NOT use GMM.
#' @param dofilter logical(1L). Filter the expression data or not.
#' @param output The output of the function. Available options: plot (default),
#' and matrix.
#' @param hclust_method The agglomeration method to be used.
#'  See \link[stats]{hclust}.
#' @param ... Parameters for \link[WGCNA]{adjacency}.
#' @return A ggplot object or a dist matrix.
#' @importFrom methods is
#' @importFrom WGCNA goodSamplesGenes adjacency TOMsimilarity
#' @export
getTOM <- function(seu,
                   softPower = 2,
                   dofilter = TRUE,
                   assay = CodexPredefined$defaultAssay,
                   output = c('plot', 'matrix'),
                   hclust_method = 'complete',
                   ...){
  stopifnot(is(seu, 'Seurat'))
  stopifnot(is.numeric(softPower))
  stopifnot(softPower>0)
  stopifnot(is.logical(dofilter))
  output <- match.arg(output)
  ## cleaning the cells
  expr <- as.data.frame(t(GetAssayData(seu,
                                       assay = assay,
                                       layer = 'data')))
  if(dofilter){
    gsg <- goodSamplesGenes(expr, verbose = 0)
    if(!gsg$allOK){
      expr <- expr[gsg$goodSamples, gsg$goodGenes, drop=FALSE]
    }
  }
  
  adjacency <- adjacency(expr, power = softPower)
  TOM <- TOMsimilarity(adjacency, verbose = 0)
  dimnames(TOM) <- dimnames(adjacency)
  if(output=='matrix'){
    return(TOM)
  }
  ggHeatmap(TOM, title='Topological overlap matrix similarity',
            hclust_method=hclust_method)
}

#' Plot Modules
#' @description
#' Plot a hierarchical clustering dendrogram and the modules by a given TOM
#' @param TOM Topological overlap matrix.
#' @importFrom stats hclust as.dist
#' @importFrom dynamicTreeCut cutreeDynamic
#' @importFrom WGCNA labels2colors plotDendroAndColors
#' @export
plotModules <- function(TOM){
  stopifnot(is.matrix(TOM))
  
  TOM.dissimilarity <- 1-TOM
  geneTree <- hclust(as.dist(TOM.dissimilarity), method = "average") 
  Modules <- cutreeDynamic(dendro = geneTree, distM = TOM.dissimilarity,
                           deepSplit = 2, pamRespectsDendro = FALSE,
                           minClusterSize = 2, verbose=0)
  names(Modules) <- geneTree$labels
  ModuleColors <- labels2colors(Modules)
  suppressWarnings(
    plotDendroAndColors(geneTree, ModuleColors, "Module",
                        dendroLabels = NULL, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05,
                        main = "Gene modules"))
  return(invisible(Modules))
}