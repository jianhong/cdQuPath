#' Get co-occurrence for given bins
#' @description
#' Calculate co-occurrence for each bin given range.
#' @param seu An seurat object
#' @param anno_col The annotation column used to define the factors.
#' @param minRadius,maxRadius Minimal and maximal radius
#' @param bins The number of the bins
#' @param spatialCoordsNames The column names of coordinates
#' @param celldist Output of \link{findNearCellsByRadius}
#' @param ... not used.
#' @return A data matrix contain the cell proportion and co-occurrence for 
#' each bin.
#' @importFrom utils combn
#' @export
#' @examples
#' # example code
#' 
coOccurrence <- function(seu,
                        anno_col = 'celltype',
                        minRadius = 0,
                        maxRadius = 1000,
                        bins = 100,
                        spatialCoordsNames = c('Cell.X', 'Cell.Y'),
                        celldist,
                        ...){
  stopifnot(is(seu, 'Seurat'))
  stopifnot(is.numeric(minRadius))
  stopifnot(is.numeric(maxRadius))
  # find the cell pairs in each bin
  breaks <- seq(minRadius, maxRadius,
                length.out = bins)
  breaks[bins+1] <- maxRadius+1
  
  # neighborhoods distance scanning
  celldist <- checkCellDistInput(seu, anno_col, celldist, maxRadius,
                                 spatialCoordsNames)
  
  celldist$bin <- findInterval(celldist$dist, breaks)
  celldist <- split(celldist[, c("cell_1", "cell_2", "anno1", "anno2")],
                    celldist$bin)
  
  cb <- getAnnoCombn(seu, anno_col)
  p <- lapply(celldist, getCelltypeCounts, 
              cellAnnoComb = cb, seu = seu, anno_col = anno_col)
  p <- lapply(p, function(.ele){
    all_p <- .ele$each$proportions ## for each
    p12 <- apply(.ele$pair$proportions_each, 1, min) ## for pair
    p1 <- vapply(cb, function(.cb) all_p[.cb[1]], FUN.VALUE = numeric(1L))
    p2 <- vapply(cb, function(.cb) all_p[.cb[2]], FUN.VALUE = numeric(1L))
    .p <- c(.ele$all$proportions, ifelse(p12==0, 0, p12/(p1*p2)))
    names(.p) <- c('cell proportion', names(cb))
    .p
  })
  p <- do.call(rbind, p)
}
