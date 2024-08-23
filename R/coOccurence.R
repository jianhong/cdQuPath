# findNearCellsByRadius <- function(
#     seu,
#     maxEuclideanDistance=sqrt(.Machine$integer.max),
#     spatialCoordsNames = c('Cell.X', 'Cell.Y')){
#   stopifnot(all(spatialCoordsNames %in% colnames(seu[[]])))
#   stopifnot(length(spatialCoordsNames)==2)
#   cell_pos_dat <- seu[[]][, spatialCoordsNames]
#   colnames(cell_pos_dat) <- c("x", "y")
#   searchcells <- RANN::nn2(data = cell_pos_dat,
#                            k = nrow(cell_pos_dat),
#                            searchtype = "radius",
#                            radius = maxEuclideanDistance)
#   names(searchcells) <- c('cell.index', 'euclideanDistance')
#   return(searchcells)
# }

getCellDistance <- function(seu, maxRadius=100,
                            spatialCoordsNames = c('Cell.X', 'Cell.Y')){
  stopifnot(is(seu, 'Seurat'))
  stopifnot(is.numeric(maxRadius))
  
  stopifnot(all(spatialCoordsNames %in% colnames(seu[[]])))
  stopifnot(length(spatialCoordsNames)==2)
  cell_pos_dat <- seu[[]][, spatialCoordsNames, drop=FALSE]
  colnames(cell_pos_dat) <- c('x', 'y')
  if(length(rownames(cell_pos_dat)!=ncol(seu))){
    rownames(cell_pos_dat) <- colnames(seu)
  }
  if(any(duplicated(rownames(cell_pos_dat)))){
    stop('duplicated cell names detected.')
  }
  ## sort the cells
  cell_pos_dat <- cell_pos_dat[order(cell_pos_dat$x, cell_pos_dat$y), ,
                               drop=FALSE]
  ## shift the scanning window in the order of x, y
  ## the window size, maxRadius
  ## output should be distance bin list with paired cells
  dist <- future_apply(cell_pos_dat, 1, function(xy, xys, maxRadius){
    pointDist <- function(x0, y0, x1, y1){
      sqrt((x1-x0)^2 + (y1-y0)^2)
    }
    xy_range <- rbind(xy-maxRadius, xy+maxRadius)
    inRangeCells <- xys[xys$x>=xy_range[1, 'x'] &
                          xys$x<=xy_range[2, 'x'] &
                          xys$y>=xy_range[1, 'y'] &
                          xys$y<=xy_range[2, 'y'], ,
                                 drop=FALSE]
    d <- pointDist(xy[1], xy[2], inRangeCells$x, inRangeCells$y)
    names(d) <- rownames(inRangeCells)
    d
  }, xys=cell_pos_dat, maxRadius=maxRadius, simplify = FALSE)
  ## remove itself
  dist <- mapply(dist, names(dist), FUN=function(.ele, .name){
    .ele[names(.ele)!=.name]
  })
  distN <- rep(names(dist), lengths(dist))
  dist <- data.frame(N1=distN,
                     N2=unlist(lapply(dist, names)),
                     distance = unlist(dist))
  rownames(dist) <- NULL
  rm(distN)
  ## remove duplicates, TODO, unit test for this
  N12 <- do.call(rbind, apply(dist[, c("N1", "N2"), drop=FALSE],
                              1, sort, simplify = FALSE))
  return(dist[!duplicated(N12), , drop=FALSE])
}

#' Get co-occurrence for given bins
#' @description
#' Calculate co-occurrence for each bin given range.
#' @param seu An seurat object
#' @param anno_col The annotation column used to define the factors.
#' @param minRadius,maxRadius Minimal and maximal radius
#' @param bins The number of the bins
#' @param spatialCoordsNames The column names of coordinates
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
                        spatialCoordsNames = c('Cell.X', 'Cell.Y')){
  stopifnot(is(seu, 'Seurat'))
  stopifnot(is.numeric(minRadius))
  stopifnot(is.numeric(maxRadius))
  stopifnot(all(spatialCoordsNames %in% colnames(seu[[]])))
  stopifnot(anno_col %in% colnames(seu[[]]))
  annos <- as.character(sort(unique(seu[[]][, anno_col])))
  if(length(annos)<2){
    stop('Less than two category in the ', anno_col, '.')
  }
  # neighborhoods distance scanning
  cd <- getCellDistance(seu=seu,
                        maxRadius=maxRadius,
                        spatialCoordsNames = c('Cell.X', 'Cell.Y'))
  # find the cell pairs in each bin
  breaks <- seq(minRadius, maxRadius,
                length.out = bins)
  breaks[bins+1] <- maxRadius+1
  cd$bin <- findInterval(cd$distance, breaks)
  cd$anno1 <- as.character(seu[[anno_col]][cd$N1, anno_col])
  cd$anno2 <- as.character(seu[[anno_col]][cd$N2, anno_col])
  cb <- combn(annos, 2, simplify = FALSE)
  names(cb) <- vapply(cb, paste, character(1L), collapse=' & ')
  cd <- split(cd[, c("N1", "N2", "anno1", "anno2")], cd$bin)
  total <- ncol(seu)
  p <- lapply(cd, function(.ele){
    # get p12, p1, p2
    all <- unique(c(.ele$N1, .ele$N2))
    all_anno <- as.character(seu[[]][all, anno_col])
    all_anno <- table(as.character(all_anno))
    all_p <- all_anno/length(all)
    each <- do.call(rbind, 
                    apply(.ele[, c("anno1", "anno2")], 
                          1, sort, simplify = FALSE))
    p12 <- vapply(cb, function(.cb){
      sum(each[, 1, drop=TRUE]==.cb[1] & each[, 2, drop=TRUE]==.cb[2])
    }, FUN.VALUE = numeric(1L))/length(all)
    p1 <- vapply(cb, function(.cb) all_p[.cb[1]], FUN.VALUE = numeric(1L))
    p2 <- vapply(cb, function(.cb) all_p[.cb[2]], FUN.VALUE = numeric(1L))
    .p <- c(length(all)/total, ifelse(p12==0, 0, p12/(p1*p2)))
    names(.p) <- c('cell proportion', names(cb))
    .p
  })
  p <- do.call(rbind, p)
}
