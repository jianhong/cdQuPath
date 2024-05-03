#' Assign cell type
#' @description
#' Assign cell types by given cell markers.
#' 
#' @param seu A Seurat object.
#' @param classifier A character vector of the cell type assignment markers. 
#' The names of the classifier is the cell type and the content of the 
#' classifier is the markers. The order of the classifier is the priority of
#' assignment for `celltype_first` column in output metadata.
#' @param method The cell type assignment method. 
#' First the function will convert the probabilities
#' to positive or negative value by cutoff 0.5.
#' Option 'Rank' will use the rank of probabilities for the positive markers. 
#' And option 'Boolean' will assign the cell type 
#' by the order of classifier for the positive markers.
#' @param cutoff The cutoff value for probabilities.
#' @param ... Not used.
#' @return A Seurat object with metadata 
#'  `celltype` (top celltype),
#'  `celltype_all` (all celltype),
#'  `celltype_first` (the first assignment by the pritority of the classifier).
#' @importFrom SeuratObject Assays Misc<-
#' @importFrom methods is
#' @export
detectCellType <- function(
    seu,
    classifier = CodexPredefined$classifier,
    method = c('Rank', 'Boolean'),
    cutoff = 0.5,
    ...
    ){
  method <- match.arg(method)
  stopifnot(is(seu, 'Seurat'))
  stopifnot('Please run fitGMM first.'=
              CodexPredefined$GMM %in% SeuratObject::Assays(seu))
  stopifnot('Please run FindClusters.'='seurat_clusters' %in% colnames(seu[[]]))
  stopifnot(is.character(classifier))
  stopifnot(length(names(classifier))==length(classifier))
  dat <- GetAssayData(seu, assay = CodexPredefined$GMM, layer = 'data')
  classifier <- classifier[classifier %in% rownames(dat)]
  stopifnot('No marker is available in Seurat object.'=length(classifier)>0)
  if(method=='Boolean'){
    dat <- dat>cutoff
  }
  dat <- dat[classifier, , drop=FALSE]
  if(method=='Boolean'){
    celltype <- apply(dat, 2, function(.e){
      names(classifier[.e])
    }, simplify = FALSE)
  }else{## method is Rank
    celltype <- apply(dat, 2, function(.e){
      w <- which(.e>cutoff)
      if(length(w)==0) w <- which(.e == max(.e, na.rm=TRUE))
      w <- names(classifier)[w]
      n <- names(classifier)[order(.e, decreasing = TRUE)]
      n[n %in% w]
    }, simplify = FALSE)
  }
  ## assign unknown cell type
  celltype[lengths(celltype)==0] <- rep('unknown', sum(lengths(celltype)==0))
  stopifnot(identical(colnames(seu), names(celltype)))
  celltype_first <- vapply(celltype, function(.e) .e[1],
                           FUN.VALUE = character(1L))
  if(method=='Boolean'){
    celltype1 <- split(celltype, seu$seurat_clusters)
  }else{
    celltype1 <- split(celltype_first, seu$seurat_clusters)
  }
  celltype1 <- lapply(celltype1, function(.e){
    N <- length(.e)
    .e <- table(unlist(.e))
    .e <- sort(.e, decreasing = TRUE)
    .e1 <- .e[names(.e) != 'unknown']
    if(length(.e1)==0){
      .e1 <- .e
    }
    N1 <- .e1[1]
    c(celltype=names(.e1)[1],
      total_count=N,
      celltype_count=N1,
      percentage=100*N1/N)
  })
  celltype1 <- do.call(rbind, celltype1)
  ct <- celltype1[, -1, drop=FALSE]
  mode(ct) <- 'numeric'
  celltype1 <- cbind(celltype=celltype1[, 1],
                     data.frame(ct))
  Misc(seu, slot='celltype') <- celltype1
  seu$celltype_full <- vapply(celltype, paste,
                              FUN.VALUE = character(1L),
                              collapse=';')
  seu$celltype_first <- celltype_first
  seu$celltype <- celltype1[seu$seurat_clusters, 'celltype']
  
  return(seu)
}
