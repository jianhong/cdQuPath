#' Do normalization 
#' @description
#' Do normalization for the Seurat object created by \link{createSeuratObj}
#' @param seu A Seurat object
#' @param normalizationMethod Normalization method for data table. Available
#' methods are 'LogNormalize', 'CLR' (centered log ratio), 'RC' (CPM), 'zscore'.
#' @param ... Parameters can be used by \link[Seurat]{NormalizeData}
#' @return A Seurat object
#' @importFrom methods is
#' @importFrom SeuratObject SetAssayData
#' @importFrom Seurat NormalizeData
#' @export
#'   
normData <- function(
    seu,
    normalizationMethod = c('LogNormalize', 'CLR', 'RC', 'zscore'),
    ...){
  stopifnot(is(seu, 'Seurat'))
  normalizationMethod <- match.arg(normalizationMethod)
  if(normalizationMethod=='zscore'){
    new.data <- scale(GetAssayData(seu,
                                   assay = CodexPredefined$defaultAssay,
                                   layer = 'counts'))
    seu <- SetAssayData(seu,
                        assay = CodexPredefined$defaultAssay,
                        layer='data',
                        new.data = new.data)
  }else{
    seu <- NormalizeData(seu, normalization.method = normalizationMethod, ...)
  }
  seu
}

#' Find variable markers
#' @description
#' Identifies markers that are used for cluster analysis by 
#' \link[Seurat]{FindVariableFeatures}
#' @param seu A Seurat object
#' @param avoidMarkers The antibody markers should be avoid. 
#' eg avoidMarkers = c('sox9', 'hist3h2a').
#' @param nFeatures Number of features to select as top variable features;
#'  only used when selection.method is set to 'dispersion' or 'vst'. See
#'  \link[Seurat]{FindVariableFeatures}.
#' @param ... other parameters used by \link[Seurat]{FindVariableFeatures}.
#' @return A Seurat object
#' @importFrom Seurat FindVariableFeatures
#' @importFrom methods is
#' @export
#' @examples
#' # example code
#' 

findVarMarkers <- function(
    seu,
    avoidMarkers,
    nFeatures = 10L,
    ...){
  stopifnot(is(seu, 'Seurat'))
  stopifnot(is.numeric(nFeatures))
  
  seu <- FindVariableFeatures(seu, ..., nfeatures = nFeatures)
  if(length(avoidMarkers)){
    meta <- seu[[CodexPredefined$defaultAssay]]@meta.data
    k <- meta$var.features %in% avoidMarkers
    status.col <- colnames(meta)[grepl("variable", colnames(meta))]
    meta[k, status.col] <- FALSE
    meta$var.features.rank[k] <- NA
    meta$var.features[k] <- NA
    seu[[CodexPredefined$defaultAssay]]@meta.data <- meta
  }
  seu
}

#' create marker correlation heatmap
#' @description
#' Create marker correlation heatmaps clustered by hclust.
#' 
#' @param seu A Seurat object
#' @param layers Layers of the default assay.
#' @param method A character string indicating which correlation coefficient
#'  is to be computed. See \link[stats]{cor}.
#' @param hclust_method The agglomeration method to be used.
#'  See \link[stats]{hclust}.
#' @param ... Parameter passed to \link[stats]{dist}
#' @return A list of ggplot object
#' @importFrom stats dist hclust cor
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot geom_tile aes scale_fill_gradientn xlab ylab 
#' theme element_text scale_x_discrete scale_y_discrete element_blank ggtitle
#' coord_fixed .data
#' @importFrom methods is
#' @importFrom SeuratObject Layers
#' @export
markerCorrelation <- function(
    seu,
    layers = Layers(seu),
    method = c('spearman', 'pearson', 'kendall'),
    hclust_method = 'complete',
    ...){
  stopifnot(is(seu, 'Seurat'))
  method <- match.arg(method)
  stopifnot(all(layers %in% Layers(seu)))
  if(length(names(layers))!=length(layers)) names(layers) <- layers
  lapply(layers, function(layer){
    tryCatch({
      mat <- GetAssayData(seu, layer = layer)
      mat <- t(as.matrix(mat))
      cor <- cor(mat, method = method)
      ggHeatmap(cor, title=layer, hclust_method=hclust_method, ...)
    }, error = function(e){
      NULL
    })
  })
}

#' Rescale data using GMM
#' @description
#' Re-scale data using multiclass Gaussian Mixture Model. If data is lack of 
#' variance, the function will trying to do log normal distribution re-scale.
#' Re-scaled values will be probabilities between 0-1.
#' @param seu A Seurat object
#' @param method The package used to do the GMM fit. Available tools are 
#' 'Rmixmod' and 'mclust'.
#' @param model If method is Rmixmod, the model will be used to fit the data.
#' see \link[Rmixmod]{mixmodGaussianModel}.
#' @param max_zero_percentage The cutoff of zero percentage. If the percentage
#' is greater than max_zero_percentage, the zeros will be removed by random
#' pick to meet the max percentage. This is used to make sure there is enough
#' variance to fit the model. 
#' @param ... Not use.
#' @return A Seurat object with a new assay named as 'GMM'. The re-scaled values
#' saved in layer 'data'. 
#' @importFrom mclust densityMclust cdfMclust
#' @importFrom Rmixmod mixmodGaussianModel mixmodCluster
#' @importFrom MASS fitdistr
#' @importFrom stats plnorm
#' @importFrom SeuratObject CreateAssayObject GetAssayData JoinLayers
#' @importFrom future.apply future_apply
#' @export
fitGMM <- function(
    seu,
    method = c('Rmixmod', 'mclust'),
    model = mixmodGaussianModel(
      family = "general",
      listModels = "Gaussian_p_Lk_Ck",
      free.proportions = FALSE, equal.proportions = TRUE
    ),
    max_zero_percentage = 0.1,
    ...){
  stopifnot(is(seu, 'Seurat'))
  method <- match.arg(method)
  dat <- GetAssayData(seu, assay = CodexPredefined$defaultAssay,
                      layer = 'counts')
  dat <- future_apply(dat, 1, function(marker){
    if(sum(marker)<2){
      # if only 1 counts or less, densityMclust will take forever.
      return(marker)
    }
    if(any(marker<0)){
      stop('Negative data detected!')
    }
    marker0 <- marker
    marker <- sort(marker)
    zero_indices <- which(marker == 0)
    zero_percentage <- length(zero_indices) / length(marker)
    if(zero_percentage>max_zero_percentage){
      toberemoved <- floor(length(marker)*
                             (zero_percentage - max_zero_percentage))
      if(toberemoved>0){
        marker <- marker[-seq.int(toberemoved)]
      }
    }
    if(method=='Rmixmod'){
      dens <- mixmodCluster(marker, nbCluster=2,
                                dataType="quantitative", models=model)
      y <- predictMixModDensity(dens, data=marker)
      y <- y[match(marker0, marker)]
    }else{
      y <- tryCatch(
        {
          dens <- densityMclust(
            marker, plot = FALSE, verbose = FALSE)
          cdf <- cdfMclust(dens, data=marker)
          cdf$y[match(marker0, marker)]
        }, error=function(.e){
          warning(.e,
                  ' Data may lack of variance.',
                  ' Trying log normal distribution rescale methods.')
          x <- marker[marker>0]
          dens <- fitdistr(x,
                           densfun = 'log-normal')
          n <- table(log(x))
          names(n) <- unique(sort(x))
          p = cumsum(n)/sum(n)
          cdf <- plnorm(p,
                        meanlog = dens$estimate['meanlog'],
                        sdlog = dens$estimate['sdlog'])
          cdf <- cdf/max(cdf)
          cdf <- c('0'=0, cdf)
          cdf[match(marker0, as.numeric(names(cdf)))]
        }
      )
    }
    y
  }, simplify = FALSE)
  dat <- do.call(rbind, dat)
  colnames(dat) <- colnames(seu)
  seu[[CodexPredefined$GMM]] <- CreateAssayObject(data = dat)
  return(seu)
}

predictMixModDensity <- function(dens, data){
  if(length(dens@results)==0){
    warning('Use simple model')
    xroot <- mean(data)
  }else{
    parameters <- dens@results[[1]]@parameters
    weight <- parameters@proportions
    mus <- as.numeric(parameters@mean)
    sigmas <- unlist(parameters@variance)
    sigmas <- sigmas[order(mus)]
    if(mus[1]==mus[2]){
      xroot <- mus[1]
    }else{
      mus <- sort(mus)
      a <- sum(c(-0.5, 0.5)/sigmas)
      b <- mus[1] / sigmas[1] - mus[2] / sigmas[2]
      c <- 0.5 * (-mus[1]^2 / sigmas[1] + mus[2]^2 / sigmas[2]) +
        log(weight[1] / weight[2]) + 0.5 * log(sigmas[2] / sigmas[1])
      xroot <- (-b - sqrt(b^2 - 4.0 * a * c)) / (2.0 * a)
    }
  }
  exp_term <- exp(data  - xroot)
  y <- exp_term / (1 + exp_term)
  return((y - min(y)) / (max(y) - min(y)))
}

#' Normalize the signal by Loess
#' @description
#' Normalize the sigal by locally polynormial regression fitting.
#' @param seu A Seurat object
#' @param markers The markers to be normalized.
#' @param pseudoCounts Numeric(1L) to avoid log2(0).
#' @param ... Parameters passed to \link{normalize.loess}.
#' @return A Seurat object with normalized values saved in layer 'counts'. 
#' The raw counts saved as a new Assay named as raw_counts.
#' @importFrom SeuratObject GetAssayData CreateAssayObject SetAssayData 
#' Embeddings
#' @importFrom Matrix sparseMatrix
#' @importFrom stats as.formula aggregate loess predict
#' @export
loessSmooth <- function(
    seu,
    markers,
    pseudoCounts = 1,
    ...){
  stopifnot(is(seu, 'Seurat'))
  if(!missing(markers)){
    stopifnot(all(markers) %in% rownames(seu))
  }else{
    markers <- rownames(seu)
  }
  dat <- GetAssayData(seu, assay = CodexPredefined$defaultAssay,
                      layer = 'counts')
  pos <- Embeddings(seu$pos)
  ## resample the positions, summary multiple points into one
  ## to reduce the data size.
  i_range <- range(pos[, 1])
  j_range <- range(pos[, 2])
  m_x <- min(1e3, length(unique(round(pos[, 1]))))
  m_y <- min(1e3, length(unique(round(pos[, 2]))))
  i <- cut(pos[, 1], breaks = m_x)
  j <- cut(pos[, 2], breaks = m_y)
  stopifnot(identical(rownames(pos), colnames(dat)))
  for(m in markers){
    message('processing marker:', m)
    mat <- aggregate(data=data.frame(value=dat[m, ],
                                     x=as.numeric(i),
                                     y=as.numeric(j)),
                     x = as.formula('value ~ x + y'),
                     FUN = function(x) log2(mean(x+pseudoCounts)))
    # test the result by ggplot(mat, aes(x, y, color=value)) + geom_point(size=1)
    mat <- sparseMatrix(i=mat$x, j=mat$y, x=mat$value,
                        dims = c(m_x, m_y),
                        dimnames=list(seq.int(m_x),
                                      seq.int(m_y)))
    mat[mat<=0] <- 1e-15 ## to avoid NaN
    # get the mean for normalize.
    mat <- normalize.loess(
      mat,
      ...)
    ## reshape the matrix to a data.frame with x, y, mean
    nmat <- mapply(i, j, FUN=function(x, y) mat[x, y])
    ndat <- dat[m, ] - nmat
    ## mask the zeros
    ndat[dat[m, ]==0] <- 0
    ## check the correlation
    # plot(as.numeric(ndat), as.numeric(dat[m, ]), type='p', cex=.2, pch=16)
    ## convert log2.it back
    dat[m, ] <- 2^ndat - pseudoCounts
    rm(mat, nmat, ndat)
  }
  ## backup the raw counts
  seu$raw_counts <- CreateAssayObject(
    counts = GetAssayData(seu, assay = CodexPredefined$defaultAssay,
                          layer = 'counts'))
  ## set the new data
  SetAssayData(seu, layer = 'counts',
               new.data = dat,
               assay = CodexPredefined$defaultAssay)
  return(seu)
}