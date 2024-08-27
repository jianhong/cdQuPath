#' Calculate Jaccard Index
#' @param seu An object of Seurat.
#' @param anno_col The annotation column used to define the factors.
#' @param maxRadius Maximal Euclidean distance of the neighbors.
#' @param spatialCoordsNames The column names of coordinates
#' @param celldist Output of \link{findNearCellsByRadius}
#' @param ... not used.
#' @return A data.frame with column names of 'A', 'B' and 'JaccardIndex'
#' @importFrom S4Vectors metadata
#' @importFrom utils combn
#' @export
JaccardIndex <- function(
    seu,
    anno_col,
    maxRadius,
    spatialCoordsNames = c('Cell.X', 'Cell.Y'),
    celldist,
    ...){
  stopifnot(is(seu, 'Seurat'))
  celldist <- checkCellDistInput(seu, anno_col, celldist, maxRadius,
                                 spatialCoordsNames)
  cmbn <- getAnnoCombn(seu, anno_col)
  cellcnt <- getCelltypeCounts(celldist, cmbn, seu, anno_col)
  ji <- lapply(cmbn, function(.ele){
    A <- cellcnt$each$counts[.ele[1]]
    B <- cellcnt$each$counts[.ele[2]]
    A_and_B <- cellcnt$pair$counts[paste(.ele, collapse = ' & '), ]
    A_and_B <- min(A_and_B)
    A_and_B/(A+B-A_and_B)
  })
  cmbn <- do.call(rbind, cmbn)
  cmbn <- data.frame(cmbn)
  colnames(cmbn) <- c('A', 'B')
  cmbn$JaccardIndex <- unlist(ji)
  return(cmbn)
}
#' Plot a heatmap for JaccardIndex output
#' @param ji output of \link{JaccardIndex}
#' @param title The title.
#' @param hclust_method The agglomeration method to be used.
#'  See \link[stats]{hclust}.
#' @param ... Parameter passed to \link[stats]{dist}
#' @return A ggplot object
#' @export
plotJaccardIndex <- function(
    ji, title = 'neighborhood Jaccard Index',
    hclust_method = 'complete', ...){
  stopifnot(is.data.frame(ji))
  stopifnot(colnames(ji)==c('A', 'B', 'JaccardIndex'))
  A <- unique(ji$A)
  cor <- matrix(0, nrow = length(A),
                ncol = length(A),
                dimnames = list(A, A))
  for(i in A){
    for(j in A){
      if(i==j){
        cor[i, j] <- 1
      }else{
        k <- ji$A==i & ji$B==j
        if(sum(k)==1){
          cor[j, i] <- cor[i, j] <- 
            as.numeric(ji$JaccardIndex[k])
        }
      }
    }
  }
  ggHeatmap(cor, title = title, hclust_method = 'complete', ...)
}