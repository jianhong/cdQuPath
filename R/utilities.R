#' check file type
#' @param path file path. character(1L)
filetype <- function(path){
  stopifnot(length(path)==1)
  if(grepl('zip$', path, ignore.case = TRUE)){
    return('unz')
  }
  f = file(path)
  on.exit(close.connection(f))
  ext = summary(f)$class
  ext
}

#' loess normalization
#' @param mat a matrix with columns containing the values of the chips
#' to normalize.
#' @param subset a subset of the data to fit a loess to.
#' @param span,family.loess parameter to be passed the function
#' \link[stats]{loess}.
#' @return the mean matrix to be subtracted.
#' @importFrom utils txtProgressBar setTxtProgressBar
normalize.loess <- function(
    mat,
    subset=sample(seq.int(nrow(mat)), min(c(5000, nrow(mat)))),
    span=2/3,
    family.loess="symmetric"){
  mat <- as.matrix(mat)
  J <- dim(mat)[2]
  II <- dim(mat)[1]
  
  w <- c(0, rep(1,length(subset)), 0) ##this way we give 0 weight to the
  ##extremes added so that we can interpolate
  means <- matrix(0,II,J) ##contains temp of what we subtract
  pb = txtProgressBar(min = 0, max = J, initial = 0, style = 3)
  for (j in 1:(J-1)){
    setTxtProgressBar(pb, j)
    for (k in (j+1):J){
      y <- mat[,j] - mat[,k]
      x <- (mat[,j] + mat[,k]) / 2
      index <- c(order(x)[1], subset, order(-x)[1])
      ##put endpoints in so we can interpolate
      xx <- x[index]
      yy <- y[index]
      aux <- NULL
      tryCatch({
        aux <- loess(yy~xx, span=span, degree=1, weights=w, family=family.loess)
        aux <- predict(aux, data.frame(xx=x), se = FALSE) / J
        means[, j] <- means[, j] + aux
        means[, k] <- means[, k] - aux
      }, warning=function(w){
        # do nothing
      }, error=function(e){
        # do nothing
      })
      rm(aux, index, x, y, xx, yy)
      #gc(reset=TRUE, verbose = FALSE)
    }
  }
  close(pb)
  ## new matrix is mat - means
  return(means)
}

#' plot heatmap by ggplot2
#' @param ht A data frame with three columns, indicates the x, y, and value
#' @param title The title
#' @param hclust_method The agglomeration method to be used.
#'  See \link[stats]{hclust}.
#' @param ... Parameter passed to \link[stats]{dist}
#' @importFrom stats dist hclust
ggHeatmap <- function(ht, title, hclust_method = 'complete', ...){
  htN <- colnames(ht)
  colv <- hclust(dist(scale(t(ht)), ...), method = hclust_method)$order
  ht <- melt(ht)
  colnames(ht) <- c('X', 'Y', 'value')
  p <- ggplot(ht, aes(x=.data$X,
                       y=.data$Y,
                       fill=.data$value)) + 
    geom_tile() +
    scale_x_discrete(limits=htN[colv]) + 
    scale_y_discrete(limits=htN[colv]) + 
    scale_fill_gradientn(colors=CodexPredefined$heatColor) +
    xlab(NULL) + ylab(NULL) + ggtitle(title) +
    coord_fixed() +
    theme(axis.text.x = element_text(angle=90, hjust=1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
  return(p)
}