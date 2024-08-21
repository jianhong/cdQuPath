#' create a KNN classifier by a given Seurat object
#' @description
#' Prepare the KNN classifier from a given Seurat object. 
#' @param seuObj An Seurat object
#' @param markerList A list of protein-marker map
#' @param doLog2transform logical(1L). Do log2(x+1) transform or not for the
#' 'RNA' assay.
#' @param assay,layer The assay and layer used to create the classifier.
#' @return An \link[e1071]{gknn} object
#' @importFrom Seurat CreateSeuratObject GetAssayData
#' @export
#' @examples
#' \dontrun{
#' markers <- readRDS(system.file('extdata', 'markers.name.map.rds',
#' package='cdQuPath'))
#' pseudo <- readRDS(system.file('extdata', 'pseudobulk.scRNAseq.rds',
#' package='cdQuPath'))
#' library(future.apply)
#' plan(multisession)
#' buildClassifierFromSeurat(pseudo, markers)
#' }
buildClassifierFromSeurat <- function(seuObj,
                                      markerList, 
                                      doLog2transform = TRUE,
                                      assay = 'RNA',
                                      layer = 'count'){
  stopifnot(is(seuObj, 'Seurat'))
  stopifnot(is.logical(doLog2transform))
  stopifnot(is.list(markerList))
  stopifnot(length(names(markerList))==length(markerList))
  ## do log2 transform
  pseudoExp <- GetAssayData(seuObj, assay = assay, layer = layer)
  if(doLog2transform[1]) pseudoExp <- log2(pseudoExp+1)
  m1 <- lapply(markerList[vapply(markerList, function(.ele) 
    any(.ele %in% rownames(pseudoExp)), FUN.VALUE = logical(1L))],
    function(.ele){
      colSums(as.matrix(pseudoExp[.ele[.ele %in% rownames(pseudoExp)], ,
                                  drop=FALSE]))
    })
  m1 <- do.call(rbind, m1)
  pseudo1 <- CreateSeuratObject(as(m1, "sparseMatrix"))
  pseudo1 <- fitGMM(pseudo1)
  gmm <- GetAssayData(pseudo1, assay = 'GMM', layer = 'data')
  gmm[is.na(gmm)] <- 0
  classifier <- buildClassifier(gmm,
                                features = rownames(gmm),
                                celltypes = colnames(gmm))
  return(classifier)
}

#' Align the single cell RNAseq data to codex data
#' @description
#' Align the scRNA-seq Seurat object to the Codex data by the classifer
#' trained from the scRNA-seq Seurat object by \link{buildClassifierFromSeurat}
#' @param codexSeu The Codex Seurat object
#' @param scSeu The single cell RAN-seq Seurat object
#' @param classifier A KNN classifer created from the scRNA-seq Seurat object by
#'  \link{buildClassifierFromSeurat}
#' @param assay,newAssayName The assay name for single cell RNAseq expressions.
#' @return An aligned Seurat object.
#' @importFrom SeuratObject Embeddings Reductions 
#'  DefaultAssay<- GetAssayData CreateAssayObject SetAssayData
#' @export
#' @examples
#' \dontrun{
#' markers <- readRDS(system.file('extdata', 'markers.name.map.rds',
#' package='cdQuPath'))
#' pseudo <- readRDS(system.file('extdata', 'pseudobulk.scRNAseq.rds',
#' package='cdQuPath'))
#' classifier <- buildClassifierFromSeurat(pseudo, markers)
#' tsv <- system.file("extdata", "test.qptiff.tsv.zip",
#'                    package = "cdQuPath",
#'                    mustWork = TRUE)
#' seu <- createSeuratObj(tsv, useValue = 'Median',
#'                        markerLocations = CodexPredefined$markerLocations)
#' seu <- fitGMM(seu)
#' alignedPseudo <- alignScRNAseqToCodex(codexSeu = seu,
#'                     scSeu = pseudo,
#'                     classifier = classifier)
#' }
alignScRNAseqToCodex <- function(codexSeu, scSeu, classifier,
                                 newAssayName='aligned',
                                 assay = 'RNA'){
  stopifnot("pos" %in% Reductions(codexSeu))
  asyN <- Layers(scSeu, assay = assay)
  if(!all(c('counts', 'data') %in% asyN)){
    stop("'counts' or 'data' must be a layer name")
  }
  codexSeu <- detectCellType(codexSeu,
                 classifier=classifier,
                 celltypeColumnName='scSeuId')
  asy <- lapply(asyN, function(i){
    x <- GetAssayData(scSeu, assay=assay, layer = i)
    x <- x[, codexSeu$scSeuId]
    colnames(x) <- colnames(codexSeu)
    x
  })
  names(asy) <- asyN
  if('count' %in% asyN){
    asy1 <- CreateAssayObject(counts=asy[['counts']])
    asy <- asy[asyN!='counts']
  }else{
    asy1 <- CreateAssayObject(counts=asy[['data']])
    asy <- asy[asyN!='data']
  }
  
  codexSeu[[newAssayName]] <- asy1
  for(i in names(asy)){
    codexSeu <- SetAssayData(codexSeu, assay = newAssayName,
                             layer = i, new.data = asy[[i]])
  }
  DefaultAssay(codexSeu) <- newAssayName
  return(codexSeu)
}