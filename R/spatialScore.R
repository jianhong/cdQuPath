#' Calculate spatial score
#' @param seu An seurat object
#' @param celltypeCenter,celltypeLeft,celltypeRight The celltype for center,
#'  left and right.
#' @param anno_col The annotation column used to define the factors.
#' @param maxRadius Minimal and maximal radius
#' @param spatialCoordsNames The column names of coordinates
#' @param celldist Output of \link{findNearCellsByRadius}
#' @param ... not used.
#' @return A data.frame with spatialScore, distance to Right celltype and
#'  distance to Left celltype.
#' @importFrom SummarizedExperiment colData
#' @importFrom S4Vectors metadata
#' @export
spatialScore <- function(
    seu,
    celltypeCenter,
    celltypeLeft,
    celltypeRight,
    anno_col = 'celltype',
    maxRadius = 1000,
    spatialCoordsNames = CodexPredefined$spatialCoordsNames,
    celldist,
    ...){
  stopifnot(is(seu, 'Seurat'))
  stopifnot(is.numeric(maxRadius))
  stopifnot(all(spatialCoordsNames %in% colnames(seu[[]])))
  
  stopifnot(is(seu, 'Seurat'))
  stopifnot(anno_col %in% colnames(seu[[]]))
  celltypes <- as.character(sort(unique(seu[[]][, anno_col])))
  stopifnot('celltypeCenter is not available'=celltypeCenter %in% celltypes)
  stopifnot('celltypeLeft is not available'=celltypeLeft %in% celltypes)
  stopifnot('celltypeRight is not available'=celltypeRight %in% celltypes)
  centerCells <- which(seu[[]][, anno_col]==celltypeCenter)
  
  celldist <- checkCellDistInput(seu, anno_col, celldist, maxRadius,
                                 spatialCoordsNames)
  
  getDist <- function(cellanno){
    dist <- celldist[celldist$anno1==celltypeCenter &
                        celldist$anno2==cellanno,
                      c('cell_1', 'dist'), drop=FALSE]
    colnames(dist)[1] <- 'cell_2'
    dist <- rbind(dist,
                   celldist[celldist$anno2==celltypeCenter &
                              celldist$anno1==cellanno,
                            c('cell_2', 'dist'), drop=FALSE])
    dist <- split(dist$dist, dist$cell_2)
    dist <- vapply(dist, min, FUN.VALUE = numeric(1L))
    dist <- data.frame(cell=as.integer(names(dist)),
                        dist=dist)
    return(dist)
  }
  distR <- getDist(celltypeRight)
  distL <- getDist(celltypeLeft)
  dist <- merge(distL, distR, by='cell', suffixes = c('.L', '.R'), all = TRUE)
  dist[is.na(dist)] <- .Machine$double.xmax
  dist$spatialScore=dist$dist.R/dist$dist.L
  dist$cell <- colnames(seu)[dist$cell]
  return(dist)
}
