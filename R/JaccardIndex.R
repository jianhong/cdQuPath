#' Calculate Jaccard Index
#' @param spe Output of \link{neighborhoodsAna}. An object of SpatialExperiment.
#' @param maxgap Maximal distance of the neighbors.
#' @param ... not used.
#' @importFrom SummarizedExperiment colData
#' @importFrom S4Vectors metadata
#' @importFrom utils combn
#' @export
JaccardIndex <- function(
    spe,
    maxgap=Inf,
    ...){
  stopifnot(is(spe, 'SpatialExperiment'))
  stopifnot(is.numeric(maxgap))
  fnc <- metadata(spe)
  stopifnot('spe must be output of neighborhoodsAna' = 
              all(c('cells', 'distance', 'anno_col', 'k') %in% names(fnc)))
  cd <- colData(spe)
  celltypes <- unique(cd[, fnc$anno_col])
  celltypes_counts <- table(cd[, fnc$anno_col])
  neighbors <- fnc$cells
  neighbors[fnc$distance>maxgap] <- NA
  cmbn <- combn(celltypes, m=2, simplify = FALSE)
  ji <- lapply(cmbn, function(.ele){
    A <- celltypes_counts[.ele[1]]
    B <- celltypes_counts[.ele[2]]
    A_and_B <- apply(neighbors[rownames(cd)[cd[,fnc$anno_col]==.ele[1]], ],
                     1, FUN = function(nc){
                       .ele[2] %in% nc
                     })
    sum(A_and_B)/(A+B-sum(A_and_B))
  })
  cmbn <- do.call(rbind, cmbn)
  cmbn <- data.frame(cmbn)
  colnames(cmbn) <- c('A', 'B')
  cmbn$JaccardIndex <- ji
  cmbn
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