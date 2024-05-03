#' Calculate spatial score
#' @param spe Output of \link{neighborhoodsAna}. An object of SpatialExperiment.
#' @param celltypeCenter,celltypeLeft,celltypeRight The celltype for center,
#'  left and right.
#' @param ... not used.
#' @return A data.frame with spatialScore, distance to Right celltype and
#'  distance to Left celltype.
#' @importFrom SummarizedExperiment colData
#' @importFrom S4Vectors metadata
#' @export
spatialScore <- function(
    spe,
    celltypeCenter,
    celltypeLeft,
    celltypeRight,
    ...){
  stopifnot(is(spe, 'SpatialExperiment'))
  stopifnot(is.character(celltypeCenter))
  stopifnot(is.character(celltypeLeft))
  stopifnot(is.character(celltypeRight))
  fnc <- metadata(spe)
  stopifnot('spe must be output of neighborhoodsAna' = 
              all(c('cells', 'distance', 'anno_col', 'k') %in% names(fnc)))
  cd <- colData(spe)
  celltypes <- unique(cd[, fnc$anno_col])
  stopifnot('celltypeCenter is not available'=celltypeCenter %in% celltypes)
  stopifnot('celltypeLeft is not available'=celltypeLeft %in% celltypes)
  stopifnot('celltypeRight is not available'=celltypeRight %in% celltypes)
  centerCells <- rownames(cd)[cd[, fnc$anno_col]==celltypeCenter]
  dist <- fnc$distance[centerCells, , drop=FALSE]
  cells <- fnc$cells[centerCells, , drop=FALSE]
  dist <- as.data.frame(t(dist))
  cells <- as.data.frame(t(cells))
  distL <- mapply(dist, cells, FUN = function(d, n){
    w <- which(n==celltypeLeft)
    if(length(w)){
      return(d[w[1]])
    }else{
      return(Inf)
    }
  })
  distR <- mapply(dist, cells, FUN = function(d, n){
    w <- which(n==celltypeRight)
    if(length(w)){
      return(d[w[1]])
    }else{
      return(Inf)
    }
  })
  return(data.frame(spatialScore=distR/distL, distR=distR, distL=distL))
}