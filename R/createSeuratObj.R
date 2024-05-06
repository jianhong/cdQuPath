#' Create Seurat object
#' @description
#' Create Seurat object from output of ExportCellDetectionMeasurement.groovy
#' @param path character(1L) or connection. The path to measurements exported by
#'  ExportCellDetectionMeasurement.groovy in QuPath.
#' @param useValue character(1L). The values should be used as count table.
#' @param markerLocations The marker locations for each marker or charcter(1L)
#' for all markers. The possible locations are 'Nucleus', and 'Cell'.
#' 'Cytoplasm' and 'Membrane' are also available but not suggested.
#' @param ... parameters could be used by \link[SeuratObject]{CreateSeuratObject}
#' @return An Seurat object.
#' @importFrom SeuratObject CreateSeuratObject CreateAssayObject 
#' CreateDimReducObject DefaultAssay Idents<-
#' @importFrom utils read.delim
#' @importFrom methods as
#' @export
#' @examples
#' # example code
#' tsv <- system.file("extdata", "test.qptiff.tsv.zip",
#'                    package = "cdQuPath",
#'                    mustWork = TRUE)
#' seu <- createSeuratObj(tsv, useValue = 'Median',
#'                        markerLocations = CodexPredefined$markerLocations)
#' seu

createSeuratObj <- function(
    path, 
    useValue = c('Mean', 'Median', 'Min', 'Max', 'Std.Dev'),
    markerLocations = 'Cell',
    ...){
  stopifnot(inherits(path, c('character', 'connection')))
  stopifnot(length(path)==1)
  stopifnot(file.exists(path))
  useValue <- match.arg(useValue)
  stopifnot(is.character(markerLocations))
  if(length(markerLocations)>1){
    stopifnot('markerLocations must use markers as vector name'=
                length(names(markerLocations))==length(markerLocations))
  }
  
  ## meta data information filter
  path_type <- filetype(path)
  if(is.character(path)){
    path <- switch(path_type,
                   'gzfile' = gzfile(path),
                   'bzfile' = bzfile(path),
                   'unz' = unz(path,
                               filename = sub('.zip$', '', basename(path),
                                              ignore.case = TRUE)),
                   'url-libcurl' = url(path),
                   path)
  }
  d <- read.delim(path)
  mcoll <- grepl(CodexPredefined$metaInfoStartWith, colnames(d))
  stopifnot('The input is not output of ExportCellDetectionMeasurement.groovy'=
              any(mcoll) && CodexPredefined$CellIDcName %in% colnames(d) &&
              CodexPredefined$multiNucleis %in% colnames(d))
  meta <- d[, mcoll, drop=FALSE]
  rownames(meta) <- meta[, CodexPredefined$CellIDcName, drop=TRUE]
  meta[[CodexPredefined$CellIDcName]] <- NULL
  meta$Cell.isMultiNucleis <- 
    ifelse(grepl('MultiNucleis', meta[[CodexPredefined$multiNucleis]]),
           'MultiNucleis', 'SingleNuclei')
  colnames(meta) <- sub('\u00b5', 'u', colnames(meta))
  colnames(meta) <- gsub('\\.+', '.', colnames(meta))
  d <- d[, !mcoll, drop=FALSE]
  cn <- colnames(d)
  genes <- sub("..Cell..Median", '', cn[grepl('Cell..Median', cn)])
  if(length(markerLocations)!=1){
    if(!all(names(markerLocations) %in% genes)){
      warning('Found markers not in the data: ',
              names(markerLocations)[!names(markerLocations) %in% genes])
      markerLocations <- markerLocations[names(markerLocations) %in% genes]
      if(length(markerLocations)<1){
        stop('No selected marker available!')
      }
    }
  }
  cn1 <- unique(sub('^(.*?)\\.\\.(Nucleus|Membrane|Cytoplasm|Cell)\\.\\.(.*?)$',
                    "..\\2..\\3", cn))
  dx <- lapply(cn1, function(.ele){
    dat <- d[, cn[grepl(.ele, cn)], drop=FALSE]
    colnames(dat) <- sub(.ele, '', colnames(dat))
    dat <- t(dat)
    colnames(dat) <- rownames(meta)
    dat
  })
  names(dx) <- sub("^_", "", gsub("\\.\\.", '_', cn1))
  dx$seurat_obj_raw <- dx[[paste0("Cell_", useValue)]]
  if(length(markerLocations)==1){
    markerLocations <- rep(markerLocations, nrow(dx$seurat_obj_raw))
    names(markerLocations) <- rownames(dx$seurat_obj_raw)
  }
  for(marker in names(markerLocations)){
    dx$seurat_obj_raw[marker, ] <- 
      dx[[paste(markerLocations[marker], useValue, sep="_")]][marker, ]
  }
  seu <- CreateSeuratObject(as(dx$seurat_obj_raw, "sparseMatrix"),
                            assay = CodexPredefined$defaultAssay,
                            meta.data = meta,
                            ...)
  dx$seurat_obj_raw <- NULL
  for(i in seq_along(dx)){
    seu[[names(dx)[i]]] <- CreateAssayObject(data = dx[[i]])
  }
  rm(dx)
  pos <- meta[, CodexPredefined$spatialCoordsNames, drop=FALSE]
  colnames(pos) <- paste('POS_', 1:2)
  seu$pos <- CreateDimReducObject(embeddings = as(pos, "matrix"),
                                  key = 'POS_',
                                  assay = DefaultAssay(seu))
  Idents(seu) <- 'Cell.isMultiNucleis'
  seu
}
