
#' @title hcrdata is a west coast abalone data-set for testing performance measures
#'
#' @description hcrdata is a fishery data-set for blacklip abalone
#'     (\emph{Haliotis rubra}) from the western zone. It constitutes five
#'     objects in a list. These are arrce a matrix of cpue from 1992 - 2020 for
#'     the 8 sau in the western zone. The yearnames, the aspirational catches
#'     from 2020. Results from a FIS = NULL, and the numbers at size = NULL.
#'
#' @name hcrdata
#'
#' @docType data
#'
#' @format A list of five objects relating to abalone fishery data
#' \describe{
#'   \item{arrce}{a 29 x 8 matrix of cpue data from 1992 - 2020}
#'   \item{yearnames}{the numbers 1992 : 2020}
#'   \item{acatches}{the aspirational catches for 2020 for sau6W to sau13W}
#'   \item{fis}{currently NULL}
#'   \item{nas}{currently NULL}
#' }
#'
#' @section Subjects:
#'  \itemize{
#'    \item performance measures
#'    \item harvest strategy
#'  }
#'
#' @source Mundy, C. and J.M. McAllister (2020) Tasmanian Abalone Assessment
#'     2019. IMAS, University of Tasmania.
#'
#' @examples
#' data(hcrdata)
#' hcrdata
NULL

