


#' @title calcexpectpopC generates population catches from TAC or acatch
#'
#' @description calcexpectpopC defines the fleet dynamics model used by the 
#'     operating model. The function converts the TAC into expected sau catches,
#'     which are then translated into population level catches while
#'     introducing management implementation error. Here this
#'     error is implemented as Log-Normal errors on diver intuitions
#'     concerning the relative abundance in each population. The error is
#'     imposed separately on the populations in each SAU. The option of TAC for
#'     a zone or aspirational catches for SAU is given to allow for different
#'     processes between jurisdictions. If only a TAC is produced by an HCR then
#'     calcexpectpopC needs to be able to generate catches by SAU prior to
#'     subdividing those SAU catches among populations. Used directly by aMSE.
#'
#' @param hcrout in Tasmania this contains the aspirational catches by SAU
#'     along with many other details.
#' @param exb the exploitable biomass at the end of the previous year. In the
#'     first year of the projections this would be the last year of the
#'     conditioning.
#' @param sauCPUE the CPUE by sau, not used by TAS but needed by SA
#' @param catsau the actual catches by sau, not used by TAS but needed by SA
#' @param fleetacatch a function that calculates catch distribution across
#'     the sau and populations. Currently not needed by Tas but needed by SA
#' @param hsargs the harvest strategy constants not used in TAS, default=hsargs
#' @param glb the globals object,  not used in TAS, default=glb 
#' @param sigmab the Log-Normal standard deviation of implementation error. The
#'     default value = 1e-08, which effectively means no errors.
#' @param year what year number are we dealing with?
#'
#' @return a vector of population catches for the year to be imposed after the
#'     estimation of the mid-year exploitable biomass (midyexpB).
#' @export
#'
#' @examples
#' print("wait on suitable internal data sets")
calcexpectpopC <- function(hcrout,exb,sauCPUE,catsau,fleetacatch,
                           hsargs=hsargs,glb=glb,sigmab=0.1,year=year) {
  # exb=zoneDP$exploitB[year-1,]; suaindex=glb$sauindex; sigmaB=0.1;
  # hcrout= list(acatch=acatch)
  sauindex <- glb$sauindex
  acatch <- hcrout$acatch
  TAC <- sum(acatch,na.rm=TRUE)  # means TAC estimation irrelevant in TAS
  totexb <- sum(exb,na.rm=TRUE)
  npop <- length(exb)
  nSAU <- length(acatch)
  error <- exp(rnorm(nSAU,mean=0,sd=sigmab)) - (sigmab^2.0)/2.0
  sauexb <- tapply(exb,sauindex,sum,na.rm=TRUE) * error
  sauexb <- sauexb * (totexb/sum(sauexb,na.rm=TRUE))
  popC <- acatch[sauindex] * (exb/sauexb[sauindex])
  popC <- popC * TAC/sum(popC)
  ans <- list(popC=popC,acatch=acatch)
  return(ans)
} # end of calcexpectpopC

#' @title catchbysau calculates the catch by sau with error for Tasmania
#'
#' @description catchbysau calculates the catch by sau with error for the
#'     Tasmanian HS. This takes the sum of the aspirational catches as the TAC,
#'     which is multiplied by the proportion of exploitable biomass for each SAU,
#'     which has log-Normal error included with sd=sigmaB.
#'
#' @param inexpB a particular years' exploitable biomass by population
#' @param sauindex a vector containing the sau index number for each population
#' @param TAC the sum of the aspirational catches for a given year
#' @param sigmab the sd of the log-normal errors included in the estimates of
#'     the exploitable biomass by SAU
#'
#' @return a vector of the catches by sau
#' @export
#'
#' @examples
#' print("wait on suitable internal data")
catchbysau <- function(inexpB,sauindex,TAC,sigmab) {
  # iter=1; year=2
  # inexpB=zoneDP$exploitB[(year - 1),,iter]
  # sauindex=sauindex;TAC=origTAC[iter]; sigmab=0.05
  totexB <- sum(inexpB,na.rm=TRUE)
  sauexpB <- tapply(inexpB,sauindex,sum,na.rm=TRUE)
  nSAU <- length(sauexpB)
  npop <- length(inexpB)
  sauexpB <- sauexpB * exp(rnorm(nSAU,mean=0,sd=sigmab))
  totsau <- sum(sauexpB,na.rm=TRUE)
  sauexpB <- sauexpB * (totexB/totsau)
  catbysau <- TAC * sauexpB/sum(sauexpB)
  return(catbysau)
} # end of catchbysau

#' @title getcpueHS recreates all the TasHS outputs from the MSE projections
#' 
#' @description getcpueHS is designed to take the cpue from the finished 
#'     projections and then for a given sau, it calculates for each year of the 
#'     projections the target cpue, the grad1score, grad4score, and targscore, 
#'     as well as the final combined score. In addition, the cpue time series 
#'     used in the calculations is also output. Used directly by aMSE.
#'
#' @param ce the (hyrs+pyrs) x nsau x reps array of cpue for the scenario found
#'     in sauout$zonePsau$cpue
#' @param catches the (hyrs+pyrs) x nsau x reps array of actual catches for the 
#'     scenario found in sauout$zonePsau$cpue 
#' @param glb the globals object
#' @param yearCE the year names for the observed CPUE used to condition the 
#'     model, found in condC$yearCE. Used to identify the observed and projected
#'     years in the simulations.
#' @param hsargs the input arguments used for the TasHS in the given scenario
#' @param sau which sau should the analysis be conducted upon
#'
#' @return a list containing the cpue and catch series used, then, for each year, the 
#'     cetarget, the grad1score, the grad4score, the target score, and the 
#'     final weighted score
#' @export
#'
#' @examples
#' print("wait on data sets")
#' # ce=cpue;catches=catches;glb=glb;yearCE=yearCE;hsargs=hsargs;sau=sau
getcpueHS <- function(ce,catches,glb,yearCE,hsargs,sau) {
  # ce=out$sauout$cpue;catches=out$sauout$catch;glb=out$glb;yearCE=out$condC$yearCE;hsargs=hsargs;sau=6
  reps <- glb$reps
  nsau <- glb$nSAU
  totyrs <- glb$hyrs + glb$pyrs
  pickyr <- match(yearCE,glb$hyrnames)
  nhyr <- length(pickyr)
  cpue <- ce[(pickyr[1]:totyrs),sau,]
  catch <- catches[pickyr[nhyr]:totyrs,sau,]
  yrs <- as.numeric(rownames(cpue))
  nyrs <- length(yrs)
  cetarg <- matrix(0,nrow=(glb$pyrs+1),ncol=reps,
                   dimnames=list(yrs[nhyr:nyrs],1:reps))
  qnt <- hsargs$targqnt
  width <- hsargs$wid
  mult <- hsargs$mult
  maxtarg <- hsargs$maxtarg[sau]
  rown <- 0
  for (targyr in nhyr:nyrs) {  # targyr=nhyr
    rown <- rown + 1
    cetarg[rown,] <- apply(cpue[1:targyr,],2,quantile,probs=qnt)
    pickT <- which(cetarg[rown,] > maxtarg)
    if (length(pickT) > 0) cetarg[rown,pickT] <- maxtarg
  }
  tmp <- apply(cpue[((nhyr-1):nyrs),],2,getgrad1)
  tmp2 <- apply(tmp,2,getscore,mult=mult)
  g1s <- tmp2[2:nrow(tmp2),]; rownames(g1s) <- yrs[nhyr:nyrs]
  tmp <- apply(cpue[((nhyr-width+1):nyrs),],2,getgrad4,wid=width)
  tmp2 <- apply(tmp,2,getscore,mult=mult)
  g4s <- tmp2[width:nrow(tmp2),]; rownames(g4s) <- yrs[nhyr:nyrs]
  if (hsargs$hcrname == "mcdahcr") {
    getts <- function(x) targscore(x,qnt=qnt,mult=mult,maxtarg=maxtarg)$scores
    targsc <- apply(cpue[(nhyr:nyrs),],2,getts); 
    rownames(targsc) <- yrs[nhyr:nyrs]
  } 
  # if (hsargs$hcrname == "constantrefhcr") {
  #   cpueyrs <- as.numeric(dimnames(ce)[[1]])
  #   if (inherits(hsargs$refperiodCE, "matrix")) { # allow for missing years
  #     refperiod <- hsargs$refperiodCE
  #   } else {
  #     refyrs <- hsargs$startCE:hsargs$endCE
  #     refperiod <- matrix(0,nrow=nsau,ncol=length(refyrs),
  #                         dimnames=list(glb$saunames,refyrs))
  #     for (i in 1:nsau) refperiod[i,] <- refyrs
  #   }    
  #   pickceyrs <- match(refperiod[sau,],cpueyrs)
  #   arrce <- ce[pickceyrs,sau,1]
  #   actualtarg <- quantile(arrce,probs=hsargs$targqnt,na.rm=TRUE)
  #   getts <- function(x) targscoreconstref(x,actualtarg,mult=mult,maxtarg=maxtarg)$scores
  #   targsc <- apply(cpue[(nhyr:nyrs),],2,getts)
  #   rownames(targsc) <- yrs[nhyr:nyrs]
  #   cetarg <- as.matrix(rep(actualtarg,(glb$pyrs+1)))
  #   rownames(cetarg) <- yrs[nhyr:nyrs]
  # }
  pmwts <- hsargs$pmwts
  finalsc <- (targsc * pmwts[1]) + (g4s * pmwts[2]) + (g1s * pmwts[3])
  index <- floor(finalsc) + 1 # add one to pick correct hcr[index]
  catchmult <- index
  for (i in 1:reps) catchmult[,i] <- hsargs$hcr[index[,i]]
  return(list(cpue=cpue,catch=catch,cetarg=cetarg,g1s=g1s,g4s=g4s,targsc=targsc,
              finalsc=finalsc,index=index,catchmult=catchmult))
} # end of getcpueHS

#' @title getgrad1 calculates the one year gradient score for all years of cpue
#'
#' @description getgrad1 calculates the one year gradient score for a
#'    vector of CPUE. The equation used is (CE(2:y) / CE(1:y-1)) - 1,
#'    which provides the annual proportional change in CPUE.
#'
#' @param vectce vector of cpue for a given spatial scale
#'
#' @return a vector of gradient1 for each year, starting from year 2
#' @export
#'
#' @examples
#'   data(hcrdata)
#'   cpue <- hcrdata$arrce
#'   nyr <- nrow(cpue)
#'   grad1 <- getgrad1(cpue[,1])
#'   score1 <- getscore(grad1)
#'   cbind(cpue[,1],grad1,score1)
getgrad1 <- function(vectce) {
  nyr <- length(vectce)
  grad1 <- c(NA,(vectce[2:nyr]/vectce[1:(nyr-1)])-1)
  grad1[is.infinite(grad1)|is.nan(grad1)] <- NA
  return(grad1)
} # end of getgrad1

#' @title getgrad4 applies a linear regression in steps of wid to input
#'
#' @description getgrad4 takes an input vector of cpue and, in chunks
#'     of length wid, converts them to proportional changes by
#'     dividing through by the first value of the short series, then
#'     applies a linear regression keeping only the gradient.
#'
#' @param vectce the input vector of cpue
#' @param wid the number of years of the input data to use in each
#'     regression
#'
#' @return a vector of length (wid-1) shorter than the input vector
#' @export
#'
#' @examples
#'   data(hcrdata)
#'   cpue <- hcrdata$arrce
#'   nyr <- nrow(cpue)
#'   x <- cpue[,1]
#'   grad4 <- getgrad4(x,wid=4)
#'   grad3 <- getgrad4(x,wid=3)
#'   cbind(x,grad4,grad3)
getgrad4 <- function(vectce,wid=4) { # vectce=ab$cpue; wid=4
  nyr <- length(vectce)
  inc <- wid-1
  num <- nyr-wid+1
  x <- 1:wid
  grad4 <- numeric(nyr)
  grad4[1:inc] <- NA
  for (i in 1:num) { # i=1
    propce <- vectce[i:(i+wid-1)]/vectce[i]
    grad4[i+inc] <- getlmcoef(propce,1:wid)[2]
    #grad4[i+inc] <- coef(lm(propce ~ x))[2]
  }
  grad4[is.infinite(grad4)|is.nan(grad4)] <- NA
  return(grad4)
} # end of getgrad4

#' @title getlmcoef is a greatly simplified way of gaining regression coefs
#'
#' @description getlmcoef is a simplified replacement for the lm function that
#'     is much faster than lm but only returns the linear regression
#'     coefficients. This is limited to a simple y ~ x model. When in use be
#'     sure to give it the values of y required.
#'
#' @param y the dependent variable
#' @param x the independent variable
#'
#' @return a vector of the intercept and gradient, in that order.
#' @export
#'
#' @examples
#' y <- c(1.226,1.237,1.238,1.156)
#' x <- 1:4
#' yp <- y/y[1]  # to make the cpue relative to the first in the series
#' getlmcoef(yp,x)
getlmcoef <- function(y,x) {
  mx <- mean(x,na.rm=TRUE)
  xres <- (x-mx)
  my <- mean(y,na.rm=TRUE)
  yres <- (y - my)
  x2 <- sum(xres^2)
  xy <- sum(xres * yres)
  grad <- xy/x2
  inter <- my - grad*mx
  return(c(inter,grad))
} # end of getlmcoef

#' @title getscore calculates the scores for the grad1 and grad4 PMs
#'
#' @description getscore calculates the scores for the grad1 and grad4
#'     performance measures. It does this by re-scaling the range of
#'     the PM values and then fitting separate linear regressions to
#'     the values above and below zero. These enable it to calculate
#'     the predicted scores.
#'
#' @param pm the raw performance measure values derived from the functions
#'     getgrad1,  getgrad4, or targscore
#' @param mult the multiplier on the bounds to expand them upwards and
#'     downwards. default value = 1.1 = 10 percent increase; from hsargs$mult
#'     
#' @seealso{
#'   \link{getgrad1}, \link{getgrad4}, \link{targscore}
#' }
#'
#' @return a vector of scores to be included in the MCDA
#' @export
#'
#' @examples
#'   data(hcrdata)
#'   cpue <- hcrdata$arrce
#'   nyr <- nrow(cpue)
#'   grad1 <- getgrad1(cpue[,1])
#'   score1 <- getscore(grad1)
#'   cbind(cpue[,1],grad1,score1)
getscore <- function(pm,mult=0.1) { # pm=grad1val[,sau]; mult=hsargs$mult
  pickb <- which((pm == (-1)) | (is.nan(pm)) | (pm == Inf) | (is.na(pm)))
  bounds <- round(extendrange(pm[-pickb],range(pm[-pickb],na.rm=TRUE),f=mult),2)
  low <- seq(bounds[1],0.0,length=6)
  high <- seq(0.0,bounds[2],length=6)
  xax <- c(low[1:5],high)
  vars <- getlmcoef(0:5,xax[1:6])
  score <- numeric(length(pm)) # just to get a vector of the correct length
  pickl0 <- which(pm <= 0)
  score[pickl0] <- pm[pickl0]*vars[2] + vars[1]
  vars2 <- getlmcoef(5:10,xax[6:11])
  pickg0 <- which(pm >= 0)
  score[pickg0] <- pm[pickg0]*vars2[2] + vars2[1]
  return(score)
} # end of getscore

#' @title makeouthcr creates an object to contain the replicate HCR outputs
#'
#' @description makeouthcr creates an objects required to hold the cumulative
#'     values for the grad4 pm. Otherwise the full run from 1992 - endyr must be
#'     recalculated each time. In this case it is a 3-D array
#'     set up to hold the values used on each SAUs previous years
#'     aspirational catch. Copies of the output array could also be used to
#'     store the scores developed for the different components of the HCR.
#'     Currently it stores 3-D arrays of years x sau x reps for g1s = gradient 1
#'     score, g4s = gradient 4 score, targsc = target CE score, finalsc = final
#'     weighted score (which can be modified by metarule 3), index = the index
#'     of the vector of hcr aspirational catch multipliers, catchmult = the 
#'     value from the hcr vector in hsargs (can be modified by metarule 1 and 2), 
#'     metaflag = a flag for when a metarule is triggered where meta rule 3 = 3,
#'     rule 2 = 2, and meta rule 1 = 1.25 (so when combined the toals can be 
#'     identifed), and cetarg
#'
#' @param glb the globals object for the MSE
#' @param hsargs the input list of arguments used by the HCR
#'
#' @return in this case a list of g1s, g4s, targsc, finalsc, index, catchmult,
#'     metaflag, and cetarg
#' @export
#'
#' @examples
#' args(makeouthcr)
makeouthcr <- function(glb,hsargs) { # glb=glb; hsargs=tashsargs
  startyr <- glb$hyrs + 1
  endyr <- glb$hyrs + glb$pyrs
  saunames <- glb$saunames
  nsau <- length(saunames)
  reps <- glb$reps
  yrnames <- c(glb$hyrnames,glb$pyrnames)
  hcryrs <- yrnames[startyr]:yrnames[endyr]
  nyrs <- length(hcryrs)
  g1s <- array(0,dim=c(nyrs,nsau,reps),
               dimnames=list(hcryrs,saunames,1:reps))
  g4s <- targsc <- finalsc <- index <- catchmult <- metaflag <- cetarg <- g1s
  outhcr <- list(g1s=g1s,g4s=g4s,targsc=targsc,finalsc=finalsc,index=index,
                 catchmult=catchmult,metaflag=metaflag,cetarg=cetarg)
  return(outhcr)
} # end of makeouthcr



#' @title targscore generates the HCR score for the target PM
#'
#' @description targscore takes in a vector of cpue x years and defines the targ
#'     from the complete series, and the limit reference point. This differs
#'     from the two getgradone and getgradwid in needing multiple years at once.
#'     One meta-rule is that the target cpue should not rise above
#'
#' @param vectce the vector of cpue to be assessed
#' @param qnt the quantile of the input vector selected as the target,
#'     default = 0.55
#' @param mult the multiplier on the bounds to expand them upwards and
#'     downwards. default value = 1.1 = 10 percent increase and decrease
#' @param maxtarg a meta-rule that sets an upper limit on the target cpue, the
#'     default=150kg/hr
#'
#' @return a list of the final year's score, the internals to the calculations,
#'     and he target and limit reference points
#' @export
#'
#' @examples
#' print("wait on suitable data")
#' args(targscore)
#' # vectce=40:140; qnt=0.55;mult=0.1; maxtarg=150.0
targscore <- function(vectce,qnt=0.55,mult=0.1,maxtarg=150.0) { #
  nyr <- length(vectce)
  targ <- quantile(vectce,probs=c(qnt),na.rm=TRUE) # adaptive ref pt.
  actualtarg <- targ  # in case targ is replaced by maxtarg
  if (targ > maxtarg) targ <- maxtarg
  bounds <- round(extendrange(vectce,range(vectce,na.rm=TRUE),f=mult),3)
  bounds[which(bounds < 0)] <- 0
  low <- seq(bounds[1],targ,length=6)
  high <- seq(targ,bounds[2],length=6)
  xax <- c(low[1:5],high)
  vars <- getlmcoef(0:5,xax[1:6])
  score <- numeric(length(vectce))
  pickl <- which(vectce <= targ)
  score[pickl] <- vectce[pickl]*vars[2] + vars[1]
  vars2 <- getlmcoef(5:10,xax[6:11])
  pickh <- which(vectce >= targ)
  score[pickh] <- vectce[pickh]*vars2[2] + vars2[1]
  rp <- c(bounds[1],targ,bounds[2],actualtarg)
  names(rp)=c("lower","target","upper","actualtarg")
  result <- tail(score,1)
  return(list(result=result,scores=score,rp=rp))
} # end of targscore

#' @title targscoreconstref generates the HCR score for the fixed target PM
#'
#' @description targscoreconstref takes in a vector of cpue x years, where the
#'     the years used are defined in hsargs as startCE and endCE, which for
#'     Tasmania default to 1992:2012 and defines the targ for each sau.
#'
#' @param vectce the vector of cpue to be assessed
#' @param realtarg the quantile of the reference period of cpue x years used as
#'     the target. The default quantile = 0.55
#' @param mult the multiplier on the bounds to expand them upwards and
#'     downwards. default value = 1.1 = 10 percent increase and decrease
#' @param maxtarg a meta-rule that sets an upper limit on the target cpue, the
#'     default=150kg/hr
#'
#' @return a list of the final year's score, the internals to the calculations,
#'     and he target and limit reference points
#' @export
#'
#' @examples
#' print("wait on suitable data")
#' args(targscore)
#' # vectce=40:140; qnt=0.55;mult=0.1; maxtarg=150.0
targscoreconstref <- function(vectce,realtarg,mult=0.1,maxtarg=150.0) { #
  nyr <- length(vectce)
  actualtarg <- realtarg
  if (realtarg > maxtarg) realtarg <- maxtarg
  bounds <- round(extendrange(vectce,range(vectce,na.rm=TRUE),f=mult),3)
  bounds[which(bounds < 0)] <- 0
  low <- seq(bounds[1],realtarg,length=6)
  high <- seq(realtarg,bounds[2],length=6)
  xax <- c(low[1:5],high)
  vars <- getlmcoef(0:5,xax[1:6])
  score <- numeric(length(vectce))
  pickl <- which(vectce <= realtarg)
  score[pickl] <- vectce[pickl]*vars[2] + vars[1]
  vars2 <- getlmcoef(5:10,xax[6:11])
  pickh <- which(vectce >= realtarg)
  score[pickh] <- vectce[pickh]*vars2[2] + vars2[1]
  rp <- c(bounds[1],realtarg,bounds[2],actualtarg)
  names(rp)=c("lower","target","upper","actualtarg")
  result <- tail(score,1)
  return(list(result=result,scores=score,rp=rp))
} # end of targscoreconstref


#' @title tasCatch extracts the catch data if required by the associated HS
#'
#' @description tasCatch no catch data are currently required by the
#'     Tasmanian HCR so this function merely returns NULL. Used directly by 
#'     aMSE.
#'
#' @param x a dummy variable not used
#'
#' @return NULL
#' @export
#'
#' @examples
#' x <- tasCatch(NA)
#' print(x)
tasCatch <- function(x) {
  return(NULL)
} # end of tasCatch

#' @title tasCPUE prepares the projected cpue data ready for the tasHCR
#'
#' @description tasCPUE prepares the projected cpue data and the yearnames from
#'    zoneDP ready for their use within the tasHCR. The yearnames identifies
#'    the extent of the years of data being included in the calculation of the
#'    cpue performance measures. In the MSE (not in reality), to ensure that 
#'    predicted values for each performance measure are available, the first 
#'    year of the predicted cpue data to be used has to be two years prior to 
#'    the first year of useful observed cpue. This is particularly necessary so 
#'    that the grad4 performance measure has four years of data to work with.
#'    Used directly by aMSE.
#'
#' @param cesau the projected cpue by SAU for a particular year and iteration
#'     from zoneDP
#' @param year the year of projection is used, along with 'startCE' to define 
#'     which years of cpue data to use.
#' @param startCE the first year of observed cpue data to use. In Tasmanian 
#'     nothing earlier than 1992 is considered sufficiently reliable for use. 
#' @param decrement how many years from the current to remove from the cpue
#'     cpue time series, default=1
#'
#' @return a list of the array of cpue and the yearnames vector
#' @export
#'
#' @examples
#' args(tasCPUE)
#' # tasCPUE(cesau=sauCPUE,year=year,startCE=hsargs$startCE,decrement=1)
#' # cesau=zoneDP$cesau; iter=1; year=59
#' # cesau=sauCPUE; year=62; startCE=hsargs$startCE; decrement=1
tasCPUE <- function(cesau,year, startCE,decrement=1) {
  yrnames <- as.numeric(rownames(cesau))
  begin <- which(yrnames == startCE)
  arrce <- cesau[begin:(year-decrement),]
  return(list(arrce=arrce,yearnames=yrnames[begin:(year-decrement)]))
} # end of tasCPUE

#' @title tasdata processes the projections to produce the data for the HCR
#'
#' @description tasdata takes the projection results from zoneDP, and histCE
#'     from 'otherdata', and appropriately samples them ready for the HCR.
#'     Used directly by aMSE. In the TasHS the no FIS data, or NaS data is used.
#'     The aCatches by SAU are used but only to be modified by the HS for the
#'     next year's acatches. tasCatch is only included because the SAHS uses
#'     catches by sau directly in HS
#'
#' @param tasCPUE a function that generates the CPUE statistics
#' @param tasFIS a function that generates the FIS statistics
#' @param tasNaS a function that generates the Numbers-at-size samples
#' @param sauCPUE the cpue by sau object from the dynamics object from the
#'     replicate model runs
#' @param sauacatch aspirational catches for each SAU from year-1, these would 
#'     be the actual catches in South Australia, so we do not need a separate
#'     tasCatch function or its equivalent in SA, as we already have the 
#'     sauacatch being input
#' @param saucatch actual catches for each SAU from year-1, these are not needed
#'     in Tasmania or Victoria, but are needed (as well as acatch) in South
#'     Australia     
#' @param sauNAS the numbers-at-size data from the zoneDP object
#' @param year the specific year being considered
#' @param decrement how many years from the current to remove from the cpue
#'     cpue time series, default=1, which means that the Aspirational catches
#'     and TAC for a year derive from the previous year (as in Tasmania). If in 
#'     SA decrement = 2 would focus on the year previous. 
#' @param startCE the initial year of available catch rates
#'
#' @return a list of arrce, yearnames, acatch, fis, and nas
#' @export
#'
#' @examples
#' print("wait on suitable data sets")
#' # sampleCE=tasCPUE;sampleFIS=tasFIS;sampleNaS=tasNAS;sauCPUE=zoneDP$cesau[,,iter]
#' # sauacatch=zoneDP$acatch[,,iter];saucatch=zoneDP$catsau[,,iter]
#' # sauNAS=list(Nt=zoneDP$Nt[,,,iter];catchN=zoneDP$catchN[,,,iter]
#' # NumNe=zoneDP$NumNe[,,,iter]);year=year;startCE=hsargs$startCE;
#' # decrement=hsargs$decrement
tasdata <- function(tasCPUE, tasFIS, tasNaS, sauCPUE, sauacatch, saucatch,
                    sauNAS, year,decrement, startCE) {
  outce <- tasCPUE(sauCPUE, year, startCE, decrement=decrement)
  yearnames <- outce$yearnames   # omit empty first year
  arrce <- outce$arrce
  acatches=sauacatch[year-1,]
  fis <- tasFIS(NA)
  nas <- tasNaS(sauNAS)
  ans <- list(arrce=arrce,yearnames=yearnames,acatches=acatches,fis=fis, nas=nas)
  return(ans)
} # end of tasdata


#' @title tasFIS calculates required FIS data
#'
#' @description tasFIS no FIS data are currently required by the Tasmanian
#'     HCR so this function merely returns NULL. Used directly by aMSE.
#'
#' @param x a dummy variable not used
#'
#' @return NULL
#' @export
#'
#' @examples
#' x <- tasFIS(NA)
#' print(x)
tasFIS <- function(x) {
  return(NULL)
}

#' @title tasNas calculates required Numbers-at-Size data
#'
#' @description tasNas no numbers-at-size data are currently required by the
#'     Tasmanian HCR so this function merely returns NULL. Used directly by 
#'     aMSE.
#'
#' @param x a dummy variable not used
#'
#' @return NULL
#' @export
#'
#' @examples
#' x <- tasNaS(NA)
#' print(x)
tasNaS <- function(x) {
  return(NULL)
}



