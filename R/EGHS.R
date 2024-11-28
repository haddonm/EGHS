

#' @importFrom utils tail
#' @importFrom stats quantile loess sd rnorm median
#' @importFrom grDevices extendrange rgb 
#' @importFrom graphics abline arrows axis points polygon grid lines barplot
#' @importFrom hplot plotprep parset getmax getmin plotnull
#' @importFrom makehtml filenametopath addplot
NULL


#' @title EGHS: functions to assist with implementing Tasmania's Abalone HS
#'
#' @description EGHS provides an array of functions that are required by the
#'     Tasmanian Abalone Harvest Strategy. These are the fixed functions. To
#'     implement the HS one still requires functions for calculating the 
#'     final Multi-Criterion Decision Analysis values with the resulting 
#'     aspirational catches by Spatial Assessment Units (which are currently
#'     the Statistical blocks). 
#'
#' @section Analysis functions:
#' \describe{
#'   \item{getgrad1}{one version for calculating the single year gradients}
#'   \item{getgrad4}{calculates the four year gradients, optional width}
#'   \item{getlmcoef}{getlmcoef is a simplified, faster way of calculating
#'                    regression coefs}
#'   \item{getscore}{calculates the scores for the grad1 and grad4 PMs}
#'   \item{targscore}{calculates the HCR score for the target PM}
#' }
#' 
#' @section The Tasmanian hsargs is a list containing:
#' \describe{
#'   \item{mult}{default=0.1, used to extend the range of observed cpue by 
#'       extending the value above and below the observed}
#'   \item{wid}{default=4; the number of years to use in grad4}
#'   \item{targqnt}{default=0.55; the quantile for the cpue target level}
#'   \item{maxtarg}{default=150; the maximum value the cpue target can reach,
#'       if the actual target lies above this it is then limited to the maxtarg}
#'   \item{acatch}{a vector of expected annual catches one value for each sau,
#'       this is only used when the consthcr hcr is used}
#'   \item{pmwts}{default=c(0.65,0.35,0.1); the relative weights given to the 
#'       three performance measures when calculating the final scores}
#'   \item{hcr}{default = C(0.25,0.75,0.8,0.85,0.9,1,1.05,1.1,1.15,1.2); the
#'       harvest control rule scales that transform the combined score into a
#'       multiplier for the previous aspirational catch for an sau given a
#'       score from 1 - 10, the score acts as an index to the vector location}
#'   \item{hcrm3}{exactly like the hcr argument but the vector of values are 
#'       those to be used if meta-rule 3 is triggered. This aims to give back
#'       more catch when the stock is improving rapidly}
#'   \item{startCE}{used in constant reference period HS, eg 2000}
#'   \item{endCE}{used in constant reference period HS, eg 2019}
#'   \item{refperiodCE}{use a vector of years instead of startCE and endCE
#'       eg 2000:2019}
#'   \item{metRunder}{default=1; if > 0 then the meta rule 1 is applied}
#'   \item{metRover}{default=1; if > 0 then the meta rule 2 is applied}
#'   \item{decrement}{default=1; working in the year after the data are 
#'       available, if decrement = 1 means use all data up to the latest year 
#'       (as in TAS). If decrement = 2 means omit the final (partial) year of 
#'       data from the assessment (as in SA}
#'   \item{pmwtSwitch}{default = 4, meta-rule 3 how many consecutive years 
#'       above the targCPUE must CPUE increase before switching the performance 
#'       measure weights from pmwts to stablewts and the hcr to hcrm3}
#'   \item{stablewts}{default= c(0.4, 0.5, 0.1); what performance measure 
#'       weights should be used once pmwtSwitch is triggered}
#'   \item{hcrname}{default="constantrefhcr" the name of the harvest control 
#'       rule used, defining a constant reference period described by 
#'       startCE and endCE, or refperiodCE. Alternives in Tasmania could be 
#'       consthcr (a constant catch HS) and mcdahcr, which implies an adaptive
#'       reference period }    
#' }
#' 
#' @name EGHS
#' @keywords internal
"_PACKAGE"
NULL
