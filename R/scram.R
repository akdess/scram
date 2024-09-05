
#' @details
#' The main functions you will need to use are CreateSCRAMObject() and runSCRAM(scram_object).
#' For additional details on running the analysis step by step, please refer to the example vignette.
#' @aliases scram-package
"_PACKAGE"


#' The scram Class
#'
#'
#' The scram Class
#'
#' @name scram
#' @rdname scram
#' @aliases scram-class
#' @exportClass scram

scram <- methods::setClass(
  "scramv2",
  slots = c(
    seurat_obj = "ANY",
    cnv = "ANY",
    snv = "ANY",
    tumor = "ANY",
    models= "ANY",
    organism= "ANY",
    min_support = "ANY",
    max_set_size = "ANY",
    cellTypeLevelAnnotation = "ANY",
    setsAnnotatedCellTypeLevel = "ANY",
     setsAnnotatedCellClassLevel = "ANY",
    cellTypeLevelAnnotationDetailed = "ANY",
    cellClassLevelAnnotation = "ANY",
    cellClassAnnotationDetailed = "ANY",
    results_list="ANY",
    casper_scell="ANY",
    casper_bulk="ANY",
    simplified_models="ANY"
  )
)
