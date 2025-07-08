#' @importFrom SpatialExperiment SpatialExperiment
#' @importFrom SingleCellExperiment SingleCellExperiment
seuToSpe <- function(
    seu,
    ...){
  stopifnot(is(seu, 'Seurat'))
  spe <- SpatialExperiment(
    assay = GetAssayData(seu, ...),
    colData = seu[[]],
    spatialCoordsNames = c('Cell.X', 'Cell.Y')
  )
}

#' Neighborhoods analysis
#' @description
#' Use hoodscanR package to do the neighborhoods analysis.
#' 
#' @param seu A Seurat object
#' @param k The maximum number of nearest cells to compute. 
#' See \link[hoodscanR]{findNearCells}.
#' @param anno_col Character vector. The name of annotation column to use.
#' See \link[hoodscanR]{findNearCells}.
#' @param clusterK The number of clusters for \link[hoodscanR]{clustByHood}
#' @param ... Other parameters used by \link[hoodscanR]{clustByHood} except
#' 'object', 'pm_cols' and 'k'.
#' @importFrom hoodscanR findNearCells scanHoods mergeHoodSpe mergeByGroup
#' calcMetrics clustByHood
#' @importFrom SeuratObject Assays
#' @importFrom S4Vectors metadata metadata<-
#' @export
neighborhoodsAna <- function(
    seu,
    k = 25,
    anno_col = 'celltype_first',
    clusterK = 12,
    ...){
  stopifnot(is(seu, 'Seurat'))
  stopifnot('Please run fitGMM first.'=
              CodexPredefined$GMM %in% SeuratObject::Assays(seu))
  spe <- seuToSpe(seu, assay = 'GMM', layer = 'data')
  # neighborhoods distance scanning, k = maximal neighbor cells
  fnc <- findNearCells(spe, k = k, anno_col = anno_col)
  # probability matrix
  pm <- scanHoods(fnc$distance)
  # Merge probability matrix based on annotations
  hoods <- mergeByGroup(pm, fnc$cells)
  # Merge probability matrix into SpatialExperiment object.
  spe <- mergeHoodSpe(spe, hoods)
  # Calculate metrics for probability matrix
  spe <- calcMetrics(spe, pm_cols = colnames(hoods))
  # Cluster the probability matrix with K-means
  spe <- clustByHood(spe, pm_cols = colnames(hoods), k=clusterK, ...)
  fnc$anno_col <- anno_col
  fnc$k <- k
  metadata(spe) <- fnc
  spe
}
