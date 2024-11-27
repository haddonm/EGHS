
#' @title checkhsargs review the hsargs and sends a summary to the console
#' 
#' @description checkhsargs reviews the arguments within hsargs and summarizes
#'     them to the console. Now that hsargs is getting relatievly complex, 
#'     mistakes are easier to make so this provides a review to ensure each
#'     scenario is being run as desired.
#'
#' @param x the hsargs object
#'
#' @return nothing 
#' @export
#'
#' @examples
#' print("wait on data")
checkhsargs <- function(x) {
  if (x$hcrname %in% c("mcdahcr","constantrefhcr","consthcr")) {
    cat("The hcr being used is: ",x$hcrname,"\n")
  } else {
    stop("Unknown hcr being used. \n")
  }
  if (x$metRunder > 0) {
    cat("meta-rule 2 is being used \n")
  }
  if (x$metRover > 0) {
    cat("meta-rule 1 is being used \n")
  } 
  if (x$pmwtSwitch > 0) {
    cat("meta-rule 3 is being used \n")
    if (sum(x$stablewts) != 1.0) cat("stablewts sum is not equal to 1! \n")
  }
  if (x$hcrname == "mcdahcr") {
    cat("An adaptive reference period is in use  \n")
  } else {
    if (inherits(x$refperiodCE,"matrix")) {
      cat("CPUE reference periods are defined as a matrix \n")
      if (length(x$endCE) > 0) 
        cat("endCE still in hsargs but is now redundant \n")
    }
    if ((length(x$endCE) > 0) & (length(x$refperiodCE) == 0)) {
      text <- paste0("hsargs$startCE and hsargs$endCE are used to define ",
                     "the reference periods for each sau \n")
      cat(text,"\n")
    }
  }  
} # end of checkhsargs

#' @title plotfinalscores plots projected catches, cpue, PMs, and final score
#'
#' @description plotfinalscores takes the individual sau output from getcpueHS
#'     and generates a plot made up of six sub-plots of the replicates for each
#'     sau. In sequence these plots are the projected catches, the projected
#'     cpue, the gradient 1 scores, the gradient 4 scores, the target cpue, and
#'     the final total score. This enables the relationships between the scores
#'     to be examined. Only used in do_MSE
#'
#' @param outhcr the HS scores objects generates by mcdahcr or constantrefhcr
#' @param zoneDP the object containing the dynamics of the projected zone
#' @param sau which sau index is to have its outputs plotted. With the 8 sau in
#'     the western zone this index is 1:8, where 1 = 6W, 2= sau7, 3=sau8, etc.
#' @param minprojC the minimum y-axis value for the projected catches
#' @param minprojCE the minimum y-axis value for the projected cpue
#' @param mintargCE the minimum y-axis value for the projected target cpue
#' @param filen default="", this would be set by the process running the 
#'     function for each sau
#'
#' @return invisibly a matrix of the median catches, cpue, targetce, finalscore,
#'     grad1 score, grad4 score, and target cpue score.
#' @export
#'
#' @examples
#' print("wait on data sets")
plotfinalscores <- function(outhcr,zoneDP,sau,minprojC=0,minprojCE=0,mintargCE=0,
                            filen="") {
  catch <- zoneDP$catsau[,sau,]
  cpue <- zoneDP$cesau[,sau,]
  yrs <- as.numeric(rownames(cpue))
  reps <- ncol(cpue)
  cetarg <- outhcr$cetarg[,sau,]
  g1s <- outhcr$g1s[,sau,]
  g4s <- outhcr$g4s[,sau,]
  targsc <- outhcr$targsc[,sau,]
  score <- outhcr$finalsc[,sau,]
  tyrs <- as.numeric(rownames(cetarg))
  pickyr <- match(tyrs,yrs)
  catch <- catch[pickyr,]
  cpue <- cpue[pickyr,]
  ymax <- getmax(catch)    # first plot catches
  medcpue <- apply(cpue,1,median)
  medcatch <- apply(catch,1,median)
  medtarg <- apply(cetarg,1,median)
  medsc <- apply(score,1,median)
  plotprep(width=8, height=9,newdev=FALSE,filename=filen,verbose=FALSE)
  parset(plots=c(7,1),margin=c(0.2,0.4,0.025,0.1))
  plot(tyrs,catch[,1],type="l",lwd=1,col="grey",ylim=c(minprojC,ymax),
       ylab="Projected Catch",xlab="",panel.first=grid(),yaxs="i")
  for (i in 1:reps) lines(tyrs,catch[,i],lwd=1,col="grey")
  lines(tyrs,medcatch,lwd=2,col=2)
  ymax <- getmax(cpue)  # max for  plot cpue
  ymin <- getmin(cpue)  # min for  plot cpue
  plot(tyrs,cpue[,1],type="l",lwd=1,col="grey",ylim=c(ymin,ymax),
       ylab="Projected CPUE",xlab="",panel.first=grid())
  for (i in 1:reps) lines(tyrs,cpue[,i],lwd=1,col="grey")
  lines(tyrs,medcpue,lwd=2,col=2)
  lines(tyrs,medtarg,lwd=2,col=4)
  # now the grad1 score
  medg1s <- apply(g1s,1,median)
  plot(tyrs,g1s[,1],type="l",lwd=1,col="grey",ylim=c(0,10),
       panel.first=grid(),ylab="Grad1 Score",xlab="")
  for (i in 1:reps) lines(tyrs,g1s[,i],lwd=1,col="grey")
  lines(tyrs,medg1s,lwd=2,col=2)
  # now the grad4 score
  medg4s <- apply(g4s,1,median)
  plot(tyrs,g4s[,1],type="l",lwd=1,col="grey",ylim=c(0,10),
       panel.first=grid(),ylab="Grad4 Score",xlab="")
  for (i in 1:reps) lines(tyrs,g4s[,i],lwd=1,col="grey")
  lines(tyrs,medg4s,lwd=2,col=2)
  # now the target cpue score
  medts <- apply(targsc,1,median)
  plot(tyrs,targsc[,1],type="l",lwd=1,col="grey",ylim=c(0,10),
       panel.first=grid(),ylab="Target CE Score",xlab="")
  for (i in 1:reps) lines(tyrs,targsc[,i],lwd=1,col="grey")
  lines(tyrs,medts,lwd=2,col=2)
  ymax <- getmax(cetarg,mult=1.01)   # now the target cpue
  ymin <- getmin(cetarg,mult=1.01)   # now the target cpue
  plot(tyrs,cetarg[,1],type="l",lwd=1,col="grey",ylim=c(ymin,ymax),
       panel.first=grid(),ylab="Target CE",xlab="")
  if (ncol(cetarg) > 1) for (i in 1:reps) lines(tyrs,cetarg[,i],lwd=1,col="grey")
  abline(h=150,lwd=2,col=2,lty=2)
  lines(tyrs,medtarg,lwd=2,col=2)
  ylabel <- "Final TasHS Score"  # now the final scores
  plot(tyrs,ceiling(score[,1]),type="l",lwd=1,col="grey",ylim=c(0,10),yaxs="i",
       panel.first=grid(),ylab=ylabel,xlab="")
  for (i in 1:reps) lines(tyrs,ceiling(score[,i]),lwd=1,col="grey")
  lines(tyrs,medsc,lwd=2,col=2)
  medscores <- cbind(medcatch,medcpue,medtarg,medsc,medg1s,medg4s,medts)
  colnames(medscores) <- c("catch","cpue","targce","finalsc","grad1","grad4",
                           "targsc")
  rownames(medscores) <- names(medcatch)
  return(invisible(medscores))
} # end of plotfinalscores

#' @title plotmultandflags plots projected catches, cpue, catchmult, and meta-flags
#'
#' @description plotmultandflags takes the individual sau output from outhcr
#'     and generates a plot made up of four sub-plots of the replicates for each
#'     sau. In sequence these plots are the projected catches, the projected
#'     cpue, the catchmult values, and a histogram of the metaflag values if 
#'     there are any meta flags raised. If no flags then a null plot is made. 
#'     Only used in do_MSE
#'
#' @param outhcr the HS scores objects generates by mcdahcr or constantrefhcr
#' @param sauans the object containing the dynamics of the projected zone 
#'     compressed to the sau level, this is sauout, which is dervied from the
#'     sauplots function
#' @param sau which sau index is to have its outputs plotted. With the 8 sau in
#'     the western zone this index is 1:8, where 1 = 6W, 2= sau7, 3=sau8, etc.
#' @param minprojC the minimum y-axis value for the projected catches
#' @param minprojCE the minimum y-axis value for the projected cpue
#' @param filen default = "", which implies the plot goes to the console.
#'
#' @return invisibly a list of the median catches, cpue, catchmult, and
#'     meta-flag
#' @export
#'
#' @examples
#' print("wait on data sets")
#' # outhcr=outhcr;sauans=sauout;sau=1;minprojC=0;minprojCE=0;filen=""
plotmultandflags <- function(outhcr,sauans,sau,minprojC=0,minprojCE=0,filen="") {
  catch <- sauans$catch[,sau,]
  cpue <- sauans$cpue[,sau,]
  depleB <- sauans$depleB[,sau,]
  yrs <- as.numeric(rownames(cpue))
  reps <- ncol(cpue)
  catchmult <- outhcr$catchmult[,sau,]
  metaflag <- outhcr$metaflag[,sau,]
  tyrs <- as.numeric(rownames(catchmult))
  pickyr <- match(tyrs,yrs)
  catch <- catch[pickyr,]
  cpue <- cpue[pickyr,]
  depleB <- depleB[pickyr,]
  ymax <- getmax(catch)    # first plot catches
  medcpue <- apply(cpue,1,median)
  medcatch <- apply(catch,1,median)
  meddepl <- apply(depleB,1,median) 
  medcm <- apply(catchmult,1,median)
  plotprep(width=8, height=8,newdev=FALSE,filename=filen,verbose=FALSE)
  parset(plots=c(5,1),margin=c(0.2,0.4,0.025,0.1))
  plot(tyrs,catch[,1],type="l",lwd=1,col="grey",ylim=c(minprojC,ymax),
       ylab="Projected Catch",xlab="",panel.first=grid(),yaxs="i")
  for (i in 1:reps) lines(tyrs,catch[,i],lwd=1,col="grey")
  lines(tyrs,medcatch,lwd=2,col=2)
  ymax <- getmax(cpue)  # max for plot cpue
  ymin <- getmin(cpue)  # min for plot cpue
  plot(tyrs,cpue[,1],type="l",lwd=1,col="grey",ylim=c(ymin,ymax),
       ylab="Projected CPUE",xlab="",panel.first=grid())
  for (i in 1:reps) lines(tyrs,cpue[,i],lwd=1,col="grey")
  lines(tyrs,medcpue,lwd=2,col=2)
  # Now the exploitable biomass depletion
  ymax <- getmax(depleB)
  plot(tyrs,depleB[,1],type="l",lwd=1,col="grey",ylim=c(0,ymax),
       ylab="Exploitable Biomass Depletion",xlab="",panel.first=grid())
  for (i in 1:reps) lines(tyrs,depleB[,i],lwd=1,col="grey")
  lines(tyrs,meddepl,lwd=2,col=2)  
  # now the catchmult score
  ymax <- getmax(catchmult)
  ymin <- getmin(catchmult)
  plot(tyrs,catchmult[,1],type="l",lwd=1,col="grey",ylim=c(ymin,ymax),
       panel.first=grid(),ylab="Catch Multiplier",xlab="")
  for (i in 1:reps) lines(tyrs,catchmult[,i],lwd=1,col="grey")
  lines(tyrs,medcm,lwd=2,col=2)
  # now the metaflag information
  if (sum(metaflag,na.rm=TRUE) > 0) {
    metaf <- as.matrix(table(metaflag))
    barplot(metaf,width=1.0,space=0.25,names.arg=rownames(metaf),beside=TRUE,
            col=c("violetred","violetred","violetred"),ylab="Frequency",
            xlab="",cex.lab=1.0)
    abline(h=0,lwd=1,col=1)
  } else {
    plotnull(msg="No meta rule flags raised")
  }
  return(invisible(list(medcatch=medcatch,medcpue=medcpue,medcm=medcm,
                        meddepl=meddepl)))
} # end of plotmultandflags


