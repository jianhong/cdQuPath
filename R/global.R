#' Predefined variables
#' @importFrom grDevices hcl.colors
#' @export
#' @examples
#' CodexPredefined$markerLocations
#' 
CodexPredefined <- list(
  spatialCoordsNames = c('Cell.X', 'Cell.Y'),
  metaInfoStartWith = '^(Nucleus|Cell)',
  CellIDcName = 'Cell.ID',
  multiNucleis = 'Cell.classification',
  defaultAssay = 'RNA',
  GMM = 'GMM',
  heatColor = hcl.colors(12, "YlOrRd", rev = TRUE),
  markerLocations = c(
    'DAPI.01'='Nucleus',
    'CD163'='Cell',
    'CD127'='Cell',
    'sox9'='Nucleus',
    'CD14'='Cell',
    'CD4'='Cell',
    'mpo'='Cell',
    'cd68'='Cell',
    'cd45'='Cell',
    'foxp3'='Cell',
    'lcp1'='Cell',
    'pd1'='Cell',
    'cd8'='Cell',
    'hif1a'='Cell',
    'ki67'='Nucleus',
    'cd19'='Cell',
    'cd206'='Cell',
    'znf90'='Cell',
    'asma'='Cell',
    'cd3'='Cell',
    'cd86'='Cell',
    'hist3h2a'='Cell',
    'cd34'='Cell',
    'CD11c'='Cell'
  ),
  classifier = c(
    "Helper T cells"             = 'CD4',
    "M1 Macrophage"              = 'cd68',
    "Cartilage tumor cell"       = 'sox9',
    "Proliferation Marker"       = 'ki67',
    "Tumor/Myeloid"              = 'lcp1',
    "MONOCYTE"                   = 'CD14',
    "M2 MACROPHAGE"              = 'cd206',
    "M2 MACROPHAGE"              = 'CD163',
    "HSC"                        = 'cd45',
    "Cytotoxic T Cells"          = 'cd8',
    "mDC/cDC - Dendritic cells"  = 'CD11c',
    "B CELLS"                    = 'cd19',
    "Endothelial"                = 'cd34',
    "Pericyte/Fibroblast"        = 'asma'
  ),
  classifier2 = list(
    positive = list(
      'B cell'                   = list( 
        symbol = c('cd45', 'cd19'), weight = c(1/2, 1/2)),
      'T cell'                   = list( 
        symbol = c('CD4'), weight = c(1)),
      'CD4_TFH'                   = list( 
        symbol = c('CD4', 'pd1'), weight = c(1/2, 1/2)),
      'T-reg'                   = list( 
        symbol = c('CD4', 'foxp3', 'pd1'), weight = c(1/3, 1/3, 1/3)),
      'CD8 T-cell'              = list( 
        symbol = c('cd8'), weight = c(1)),
      'Dendritic cell'          = list( 
        symbol = c('CD11c'), weight = c(1)),
      'NK'                      = list( 
        symbol = c('cd57'), weight = c(1)),
      'Macrophage'              = list( 
        symbol = c('cd68'), weight = c(1)),
      'M2 macrophage'           = list( 
        symbol = c('cd68', 'CD163'), weight = c(1/2, 1/2)),
      'HSC'                     = list( 
        symbol = c('cd45'), weight = c(1)),
      'cd45+lcp1'               = list( 
        symbol = c('cd45', 'lcp1'), weight = c(1/2, 1/2)),
      'Monocyte'                = list( 
        symbol = c('CD14'), weight = c(1)),
      'Stromal'                = list( 
        symbol = c('cd34', 'asma'), weight = c(1/2, 1/2)),
      'Tumour'                = list( 
        symbol = c('sox9', 'ki67'), weight = c(1/2, 1/2))
    ),
    negative = list(
      'Tumour'                = list(
        symbol = c('cd45'), weight = c(1))
    )
  )
)

