buildFeatureVector <- function(classifier, data){
  cn <- attr(classifier$term, 'term.labels')
  ## rownames of data 
  if(all(cn %in% rownames(data))){
    data <- t(data)
  }else{
    if(!all(cn %in% colnames(data))){
      stop("classifier is not match the features in the data.")
    }
  }
  data[, cn]
}

#' Build a generalized k-Nearest Neighbors Classifier for given reference
#' @description
#' For a given features data, comutes the k-nearest neighbours classifier
#' making use of general distance.
#' @param data A data.frame, containing features for \link[e1071]{gknn}
#'  training. Data must be normalized into the range of 0 to 1.
#' @param features The features used for cell type prediction.
#' @param celltypes The cell types for each row or column of the data.
#' @return An \link[e1071]{gknn} model.
#' @importFrom e1071 gknn
#' @export
#' @examples
#' data(iris)
#' iris[, 1:4] <- apply(iris[, 1:4], 2, function(.ele){
#'   (.ele-min(.ele))/max(.ele)
#' })
#' classifier <- buildClassifier(iris, colnames(iris)[1:4], iris$Species)
#' predictCelltypes(classifier, iris[c(1, 51, 101), ])
buildClassifier <- function(data, features, celltypes){
  if(all(features %in% rownames(data))){
    data <- t(data)
  }else{
    if(!all(features %in% colnames(data))){
      stop("Not all features in the input data")
    }
  }
  if(length(celltypes)!=nrow(data)){
    stop('The length of celltypes do not equal to the length of the data.')
  }
  if(max(data[, features])>1.01 || min(data[, features])< -0.01){
    stop('The data must be normalized to [0, 1]')
  }
  data <- data.frame(data[, features], celltype=celltypes)
  model <- gknn(celltype ~., data = data)
  return(model)
}

#' predict celltypes by gknn classifier
#' @description classify putative cell types by given classifier.
#' @param classifier a gknn classifier exported by \link{buildClassifier}.
#' @param data A data.frame, containing features for \link[e1071]{gknn}
#'  prediction.
#' @return Celltype factors for each cell.
#' @export
#' @importFrom stats predict
#' @examples
#' data(iris)
#' iris[, 1:4] <- apply(iris[, 1:4], 2, function(.ele){
#'   (.ele-min(.ele))/max(.ele)
#' })
#' classifier <- buildClassifier(iris, colnames(iris)[1:4], iris$Species)
#' predictCelltypes(classifier, iris[c(1, 51, 101), ])
#' 
predictCelltypes <- function(classifier, data){
  data <- buildFeatureVector(classifier = classifier,
                             data = data)
  celltypes <- predict(classifier, data, type="class")
  return(celltypes)
}