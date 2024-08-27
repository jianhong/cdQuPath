#' Find the nearby cells by a given distance
#' @description
#' Using RANN::nn2 function to get the cell-cell distance by a given distance.
#' @param seu An seurat object.
#' @param maxRadius numeric(1L). Maximal Euclidean distance.
#' @param spatialCoordsNames The column names of coordinates
#' @return A data.frame of distance pairs with column name
#'  as 'cell_1', 'cell_2', 'dist'
#' @importFrom RANN nn2
#' @export
#' @examples
#' # example code
#' 
findNearCellsByRadius <- function(
    seu,
    maxRadius,
    spatialCoordsNames = c('Cell.X', 'Cell.Y')){
  stopifnot(is(seu, 'Seurat'))
  stopifnot(is.numeric(maxRadius))
  stopifnot(length(maxRadius)==1)
  stopifnot(all(spatialCoordsNames %in% colnames(seu[[]])))
  stopifnot(length(spatialCoordsNames)==2)
  cell_pos_dat <- seu[[]][, spatialCoordsNames]
  colnames(cell_pos_dat) <- c("x", "y")
  res <- RANN::nn2(data = cell_pos_dat,
                   k = nrow(cell_pos_dat),
                   searchtype = "radius",
                   radius = maxRadius)
  names(res) <- c('cell.index', 'euclideanDistance')
  res <- getCellDistance(res)
  return(res)
}

#' @noRd
#' @param res output of \link{findNearCellsByRadius}
#' @param undefinedIdx The default index number if not find. 
#' (see \link[RANN]{nn2})
#' @param undefinedDist The default distance if not find.
#' (see \link[RANN]{nn2})
getCellDistance <- function(res, undefinedIdx=0, undefinedDist=1.340781e+154){
  stopifnot(names(res)==c('cell.index', 'euclideanDistance'))
  stopifnot(identical(dim(res[['cell.index']]),
                      dim(res[['euclideanDistance']])))
  pairsIdx <- seq.int(nrow(res[['cell.index']]))
  keep <- !duplicated(res[['cell.index']])
  res <- lapply(res, function(.ele) .ele[keep, , drop=FALSE])
  pairsIdx <- pairsIdx[keep]
  pairsIdx <- rep(pairsIdx, each=ncol(res[['cell.index']]))
  pairs <- data.frame(
    cell_1 = pairsIdx,
    cell_2 = as.integer(t(res[['cell.index']])),
    dist = as.numeric(t(res[['euclideanDistance']]))
  )
  pairs <- pairs[pairs[, 'cell_2']!=undefinedIdx, , drop=FALSE]
  pairs <- pairs[pairs[, 'cell_1']!=pairs[, 'cell_2'], , drop=FALSE]
  pairs <- unique(pairs)
  ## remove duplicates again
  pairs[, c(1, 2)] <- do.call(rbind,
                              apply(pairs[, c(1, 2)],
                                    1, sort, simplify = FALSE))
  pairs <- unique(pairs)
  return(pairs)
}

#' @noRd
getAnno <- function(seu, anno_col, celldist){
  stopifnot(is(seu, 'Seurat'))
  stopifnot(anno_col %in% colnames(seu[[]]))
  stopifnot(is.data.frame(celldist))
  celldist$anno1 <- as.character(seu[[anno_col]][celldist$cell_1, anno_col])
  celldist$anno2 <- as.character(seu[[anno_col]][celldist$cell_2, anno_col])
  return(celldist)
}

#' @noRd
#' @importFrom utils combn
getAnnoCombn <- function(seu, anno_col){
  stopifnot(is(seu, 'Seurat'))
  stopifnot(anno_col %in% colnames(seu[[]]))
  annos <- as.character(sort(unique(seu[[]][, anno_col])))
  if(length(annos)<2){
    stop('Less than two category in the ', anno_col, '.')
  }
  cellAnnoComb <- combn(annos, 2, simplify = FALSE)
  names(cellAnnoComb) <- vapply(cellAnnoComb, paste,
                                character(1L), collapse=' & ')
  return(cellAnnoComb)
}

#' @noRd
getCelltypeCounts <- function(celldist, cellAnnoComb, seu, anno_col){
  # get p12, p1, p2
  stopifnot(is.data.frame(celldist))
  stopifnot(all(c('cell_1', 'cell_2', 'anno1', 'anno2') %in%
                  colnames(celldist)))
  all <- unique(c(celldist$cell_1, celldist$cell_2))
  all_anno <- as.character(seu[[]][all, anno_col])
  names(all_anno) <- as.character(all)
  all_anno_cnt <- table(as.character(all_anno))
  all_p <- all_anno_cnt/length(all)
  count12 <- split(celldist[, c('cell_1', 'cell_2')],
                   celldist[, c('anno1', 'anno2')],
                   sep = ' & ')
  count12 <- lapply(count12, unlist, use.names = FALSE)
  count12 <- lapply(cellAnnoComb, function(.ele){
    cellIdx <- unique(c(count12[[paste(.ele, collapse = ' & ')]],
                     count12[[paste(rev(.ele), collapse = ' & ')]]))
    table(all_anno[as.character(cellIdx)])[.ele]
  })
  p1_2 <- unlist(lapply(count12, function(.ele){
    if(length(names(.ele))==0) return(0)
    max(c(0, min(.ele/all_anno_cnt[names(.ele)], na.rm=TRUE)), na.rm=TRUE)
  }))
  count12 <- do.call(rbind, count12)
  colnames(count12) <- c('anno1', 'anno2')
  p12 <- count12/length(all)
  p12[is.na(p12)] <- 0
  ## list of list
  res <- list(
    'all' = list(
      counts = length(all),
      proportions = length(all)/ncol(seu)
    ),
    'each' = list(
      counts = all_anno_cnt,
      proportions = all_p
    ),
    'pair' = list(
      counts = count12,
      proportions = p1_2,
      proportions_each = p12
    )
  )
  return(res)
}

checkCellDistInput <- function(seu, anno_col, celldist, maxRadius,
                               spatialCoordsNames){
  stopifnot(is.numeric(maxRadius))
  stopifnot(all(spatialCoordsNames %in% colnames(seu[[]])))
  if(missing(celldist)){
    celldist <- findNearCellsByRadius(
      seu=seu,
      maxRadius=maxRadius,
      spatialCoordsNames = spatialCoordsNames)
  }else{
    stopifnot('celldist must be the output of findNearCellsByRadius'=
                all(c("cell_1", "cell_2", "dist") %in% colnames(celldist)) &&
                is.data.frame(celldist))
  }
  celldist <- getAnno(seu, anno_col, celldist)
  return(celldist)
}
