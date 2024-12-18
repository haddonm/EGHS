

#' @title constantrefhcr conducts the MCDA using a constant reference period
#'
#' @description constantrefhcr generates a cpue target reference point using a
#'     fixed reference period, which is input using hsargs. Just like mcdahcr it
#'     conducts the MCDA on the basis of a vector of cpue and
#'     other details prescribed in the function's arguments. Most importantly,
#'     it returns the aspirational catches by sau, the TAC for the zone, and 
#'     there are other details. The final score is rounded up to the nearest 
#'     integer and that is the index within the harvest control rule.
#'     Thus a score of <=1 points to the first cell, >1 and <=2 points to the
#'     second cell, and so on, up to a score between >9 and <=10 which points
#'     to the last cell. All the arguments used within the mcda are brought in,
#'     inside a list named hsargs, and hence are flexible. This list includes
#'     the 'wid', the number of years to use in grad4, default=4; the
#'     'targqnt', the quantile for the cpue target level, default=0.55, the
#'     'pmwts', performance measure weights, default = 0.65, 0.25, 0.1 for the
#'     target, grad4, and grad1 PMs respectively, so the order matters; the
#'     'hcr', the harvest control rule scales that transform the combined score
#'     into a TAC multiplier. A vector of 1 - 10 where each cell index
#'     represents the upper limit of the the combined score. Now the meta rules 
#'     have been included, set them to zero to exclude them. Slightly modified 
#'     to allow for single SAU.
#'
#' @param indat a list of the matrix of cpue x years that forms the
#'     basis of the assessment, the yearnames for the years making up
#'     the cpue data, the aspirational catches from the previous year (that
#'     will be multiplied to give the coming year's aspirational catches), and
#'     a NULL fis and NaS value
#' @param hsargs the arguments used within the Harvest strategy's HCR. See the
#'     description for details.
#' @param glb the globals object, this needs to contain glb$saunames and 
#'     glb$hyrs the number of years of cpuie data
#' @param projyear an index of the year for which the acatch are being estimated
#'     thus, if one has 29 years of data hen projyear would be 30.
#' @param outhcr a list of arrays generated by makehcrout or makeouthcr from the
#'     TasHS package, depending on whether you are in the projection or do_MSE 
#'     functions, within the aMSE package. It becomes the repository of all 
#'     the HS scores and other outputs.
#' @param iter which iteration within reps is being run. Needed to output the 
#'     scores to outhcr
#'     
#' @seealso {
#'     \link{getgrad1}, \link{getgrad4}, \link{getscore},
#'     \link{targscoreconstref}, \link{getlmcoef}
#' }  
#'
#' @return a list of the acatch, TAC, TAC multiplier, the score, all the
#'     details and a matrix of the reference points
#' @export
#'
#' @examples
#'   data(hcrdata)
#'   iter <- 1
#'   hsargs <- list(mult=0.1,wid = 4,targqnt = 0.55,
#'             maxtarg = c(150,150,150,150,150,150,150,150),
#'             pmwts = c(0.65,0.25,0.1),
#'             hcr = c(0.25,0.75,0.8,0.85,0.9,1,1.05,1.1,1.15,1.2),
#'             hcrm3 = c(0.25,0.75,0.8,0.85,0.9,1,1.1,1.15,1.2,1.25),
#'             startCE = 1992,endCE = 2020,metRunder = 2,metRover = 2,
#'             decrement=1, pmwtSwitch = 4,stablewts = c(0.4, 0.5, 0.1),
#'             hcrname="mcdahcr")
#'     glb <- list(saunames=c("sau6","sau7","sau8","sau9","sau10","sau11",
#'                            "sau12","sau13"),yrnames=hcrdata$yearnames,
#'                 hyrs=29,pyrs=1,reps=1,hyrnames=1992:2020,pyrnames=2021)
#'   outhcr <- makeouthcr(glb,hsargs)          
#'   hcrout <- constantrefhcr(hcrdata,hsargs=hsargs,glb,
#'                            projyear=(nrow(hcrdata$arrce)+1),outhcr,iter)
#'   str(hcrout)
constantrefhcr <- function(indat,hsargs,glb,projyear,outhcr,iter) {
  # indat=hcrdata;hsargs=hsargs;glb=glb;projyear=(nrow(hcrdata$arrce)+1);outhcr=outhcr;iter=iter  
  saunames <- glb$saunames
  arrce <- as.matrix(indat$arrce)
  nsau <- ncol(arrce)
  yrce <- nrow(arrce)
  pmwts <- hsargs$pmwts
  stablewts <- hsargs$stablewts
  pmwtSwitch <- hsargs$pmwtSwitch
  metRunder <- hsargs$metRunder
  metRover <- hsargs$metRover
  yearnames <- indat$yearnames
  acatches <- indat$acatches
  # define storage matrices
  grad1val <- matrix(0,nrow=yrce,ncol=nsau,dimnames=list(yearnames,saunames))
  grad4val <- targval <- score1 <- score4 <- scoret <- scoretot <- grad1val
  multTAC <- indexI <- grad1val
  refpts <- matrix(0,nrow=nsau,ncol=4,
                   dimnames=list(saunames,c("low","trp","high","realtrp")))
  if (length(hsargs$refperiodCE) > 0) { # allow for missing years
    refperiod <- hsargs$refperiodCE
  } else {
    refyrs <- hsargs$startCE:hsargs$endCE
    refperiod <- matrix(0,nrow=nsau,ncol=length(refyrs),
                        dimnames=list(saunames,refyrs))
    for (i in 1:nsau) refperiod[i,] <- refyrs 
  }
  flagmeta <- matrix(0,nrow=yrce,ncol=nsau,dimnames=list(yearnames,saunames))  
  for (sau in 1:nsau) {  #  sau=1
    pickceyrs <- match(refperiod[sau,],indat$yearnames)
    actualtarg <- apply(as.matrix(arrce[pickceyrs,]),2,quantile,probs=hsargs$targqnt,
                        na.rm=TRUE)
    pickce <- which(!is.na(arrce[,sau]))
    tmp <- getgrad1(arrce[pickce,sau])                     # grad1
    nec <- length(tmp)
    if (nec < yrce) tmp <- c(rep(NA,(yrce-nec)),tmp)
    grad1val[,sau] <- tmp
    score1[,sau] <- getscore(grad1val[,sau],mult=hsargs$mult)
    tmp2 <- getgrad4(arrce[pickce,sau],wid=hsargs$wid)      # grad4
    nec2 <- length(tmp2)
    if (nec2 < yrce) tmp2 <- c(rep(NA,(yrce-nec2)),tmp2) # allow for 6 and 13
    grad4val[,sau] <- tmp2
    score4[,sau] <- getscore(grad4val[,sau],mult=hsargs$mult)
    # Now estimate target from fixed period using targscoreconstref
    tmp3 <- targscoreconstref(arrce[pickce,sau],actualtarg[sau],
                              mult=hsargs$mult,maxtarg=hsargs$maxtarg[sau])
    nec3 <- length(pickce)
    scrs <- tmp3$scores
    if (nec3 < yrce) scrs <-  c(rep(NA,(yrce-nec3)),scrs)
    scoret[,sau] <- scrs
    targval[pickce,sau] <- arrce[pickce,sau]
    scoretot[, sau] <- pmwts[1] * scoret[, sau] + pmwts[2] * score4[, sau] +
      pmwts[3] * score1[, sau]
    ## Meta Rules section -----  
    mr3fire <- FALSE
    if (pmwtSwitch > 0) { # meta rule 3  pmwt-----
      # Metarule to switch weights & update scoretot
      if (all(tail(arrce[pickce, sau], pmwtSwitch) > refpts[sau, 2])) {
        tmpmr3 <- diff(arrce[, sau])
        if (all(tail(tmpmr3, pmwtSwitch) > 0)) {
          mr3fire <- TRUE
          scoretot[, sau] <- stablewts[1] * scoret[, sau] + 
            stablewts[2] * score4[, sau] + stablewts[3] * score1[, sau]
          flagmeta[yrce,sau] <- 3
        }
      } 
    } # meta rule 3 part 1 finished------
    pickI <- floor(scoretot[,sau]) + 1 # add one to get correct hcr[index]
    pickI[pickI < 1] <- 1
    pickI[pickI > 10] <- 10
    if (pmwtSwitch > 0) {   # meta rule 3 - part 2  hcr shift CM: add----
      multTAC[, sau] <- hsargs$hcrm3[pickI]
    } else {
      multTAC[, sau] <- hsargs$hcr[pickI]
    } # meta rule 3 part 2 hcr finished-----
    refpts[sau,] <- tmp3$rp # the only reference point is the target cpue
    # Meta Rule 1 starts  overs-----
    if ((metRover > 0) & (flagmeta[yrce,sau] != 3)) {     
      if (all(tail(arrce[, sau], metRover) > refpts[sau, 2])) { #CE > Target for 2 years
        tmp5 <- diff(arrce[, sau])
        if (all(tail(tmp5, metRover) > 0)) {
          multTAC[yrce, sau] <- hsargs$hcr[pickI[yrce]] # already in place
        } else {
          multTAC[yrce, sau] <- 1 # if above but not increasing assign 1 to multTAC for no change
          flagmeta[yrce,sau] <- 1
        }
      }
    }  # end of meta rule 1--------------
    if (metRunder > 0) {  # Meta Rule 2 starts under--------------
      if (all(tail(arrce[, sau], metRunder) < refpts[sau, 2])) {
        tmp6 <- diff(arrce[, sau])
        if (all(tail(tmp6, metRunder) >= 0))  {
          multTAC[yrce, sau] <-  1 # if below and increasing 2 years assign 1 to multTAC for no change
          flagmeta[yrce, sau] <- 2
        } else {
          multTAC[yrce, sau] <- hsargs$hcr[pickI[yrce]] # already in place
        }
      }
    } # end meta rule 2---------------
    indexI[,sau] <- pickI
  } # end of sau loop
  startyr <- glb$hyrs + 1
  inyr <- projyear - startyr + 1
  outhcr$g1s[inyr,,iter] <- score1[yrce,]
  outhcr$g4s[inyr,,iter] <- score4[yrce,]
  outhcr$targsc[inyr,,iter] <- scoret[yrce,]
  outhcr$finalsc[inyr,,iter] <- scoretot[yrce,]
  outhcr$index[inyr,,iter] <- indexI[yrce,]
  outhcr$catchmult[inyr,,iter] <- multTAC[yrce,]
  outhcr$metaflag[inyr,,iter] <- flagmeta[yrce,]
  outhcr$cetarg[inyr,,iter] <- refpts[,"trp"]
  acatch <- acatches * multTAC[yrce,] # to give whole numbers
  TAC <- sum(acatch,na.rm=TRUE)
  details <- list(grad4=grad4val,grad1=grad1val,targval=targval)
  out <- list(acatch=acatch,TAC=TAC,outhcr=outhcr,details=details,refpts=refpts)
  return(out)
} # end of constantrefhcr
