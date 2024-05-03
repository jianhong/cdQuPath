#' Calculate the dynamic range.
#' @description
#' The standard approach assessing dynamic range is to calculate a 
#' signal-to-background (SNR) ratio by dividing the average of the top 20 
#' brightest cells by the average intensity of the weakest 10% of cells.
#' An SNR of 10 or more supports reliable image analysis.
#' The recommend SNR range is (10, Inf], typically > 100.
#' Lower than 3 indicates poor-performing antibodies.
#' @param seu A Seurat object.
#' @param topN,bottomN A numeric.
#'  If it is greater than 1, use top topN brightest cells or bottomN weakest
#'  cells.
#'  If it is smaller than 1 and greater than 0, use top 100*topN/bottomN
#'   percentage cells.
#' @param layer The matrix used to calculate the dynamic range.
#' @param ... parameters for GetAssayData
#' @return 
#' A list with element brightest, weakest, and dynamic_range.
#' For dynamic ranges, Inf indicates zero background. NaN indicates zero/zero
#' @importFrom SeuratObject GetAssayData
#' @importFrom methods is
#' @export
#' @examples
#' # example code
#' 

dynamicRanges <- function(
    seu,
    topN = 20,
    bottomN = 0.1,
    layer = 'counts',
    ...){
  stopifnot(is(seu, 'Seurat'))
  N <- ncol(seu)
  stopifnot(is.numeric(topN))
  stopifnot(is.numeric(bottomN))
  if(topN>N || topN<=0){
    stop('topN is greater than available cells or smaller than 0.')
  }
  if(bottomN>N || bottomN<=0){
    stop('bottomN is greater than available cells or smaller than 0.')
  }
  if(topN>=1){
    topN <- seq.int(topN)
  }else{
    if(topN>0){
      topN <- seq.int(ceiling(N*topN))
    }
  }
  if(bottomN>=1){
    bottomN <- seq(N-bottomN, N)
  }else{
    if(bottomN>0){
      bottomN <- seq.int(N)[-seq.int(floor(N*bottomN))]
    }
  }
  dx <- GetAssayData(seu, layer = layer, ...)
  t20 <- 
    apply(dx, 1, function(x){
      x <- sort(x, decreasing = TRUE)
      mean(x[topN], na.rm=TRUE)
    }, simplify = TRUE)
  b10 <- 
    apply(dx, 1, function(x){
    x <- sort(x, decreasing = TRUE)
    mean(x[bottomN], na.rm=TRUE)
  }, simplify = TRUE)
  dynamicRange = t20/b10
  list(brightest = t20,
       weakest = b10,
       dynamic_range = dynamicRange)
}
