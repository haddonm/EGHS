



#' @title ploteHSphase generates a phase plot of CPUE Target vs Gradient scores
#' 
#' @description The empirical Tasmanian harvest strategy uses different CPUE
#'     statistics to produce scores that represent different facets of CPUE.
#'     Plotting the score for the CPUE gradient measure against the CPUE target 
#'     measure is used as a proxy for a phase plot of fishing mortality against 
#'     mature biomass. Strictly, statistics relating to CPUE are concerned with
#'     the exploitable biomass, which, because of the LML, cam be markedly
#'     different from the spawning biomass. Nevertheless, trend in one tend to
#'     be similar to trends in the other.
#'
#' @param cpue a vector of cpue data either predicted or observed.
#' @param rundir the directory for all files and results, default=""
#' @param wid default=4 the number of years used on getgrad4 from TasHS
#' @param fnt default = 7 = times bold. 
#' @param cex default font size = 1.0
#' @param console should the plot go to the console or be saved? Default=TRUE  
#'
#' @return nothing but it does plot a graph
#' @export
#'
#' @examples
#' print("wait on a data set")
ploteHSphase <- function(cpue,rundir="",wid=4,fnt=7,cex=1.0,console=TRUE) {
  g4 <- getgrad4(cpue)
  sc4 <- getscore(g4)
  g1 <- getgrad1(cpue)
  sc1 <- getscore(g1)
  ts <- targscore(cpue)$scores
  pickp <- which(sc4 > 0)
  years <- as.numeric(names(cpue[pickp]))
  cedat <- cbind(sc1[pickp],sc4[pickp],ts[pickp])
  rownames(cedat) <- years
  colnames(cedat) <- c("grad1",'grad4','cetarg')
  n <- nrow(cedat)
  ts <- cedat[,"cetarg"]
  sc4 <- cedat[,"grad4"]-5
  yrge <- range(years)
  if (console) { filen <- "" } else {
    filen <- filenametopath(rundir,"eHS_phase_summary.png")
  }
  plotprep(width=7,height=4.5,newdev=FALSE,filename=filen,verbose=FALSE,
           usefont=fnt)
  parset(margin = c(0.45, 0.45, 0.1, 0.1),cex=cex,cex.lab=(cex*1.1))
  plot(ts,sc4,type="p",pch=16,xlim=c(0,10),ylim=c(5,-5),xaxs="i",yaxs="i",
       xlab=paste0("eHS CPUE Target Score - Abundance Proxy  ",
                   yrge[1]," - ",yrge[2]), xaxt="n",yaxt="n",
       ylab="eHS CPUE Gradient Score - Fishing Mortality Proxy")
  axis(side=1,at=0:10)
  axis(side=2,at=seq(5,-5,-1))
  polygon(x=c(0,1,1,0,0),y=c(0,0,-5,-5,0),col=rgb(255,0,0,175,maxColorValue=255))
  polygon(x=c(1,10,10,1,1),y=c(0,0,-5,-5,0),col=rgb(255,255,0,170,maxColorValue=255))
  polygon(x=c(0,1,1,0,0),y=c(5,5,0,0,5),col=rgb(255,255,0,175,maxColorValue=255))
  polygon(x=c(1,10,10,1,1),y=c(5,5,0,0,5),col=rgb(0,255,0,120,maxColorValue=255))
  abline(v=c(1,5),lwd=2,lty=c(1,2))
  abline(h=0,lwd=2)
  points(ts,sc4,pch=16,cex=1.5)
  arrows(x0=ts[1:(n-1)],y0=sc4[1:(n-1)],x1=ts[2:n],y1=sc4[2:n],lwd=2,length=0.15)
  if (!console) {
    caption <- "eHS Phase plot of CPUE Gradient vs CPUE Target."
    addplot(filen,rundir,category="fishery",caption=caption)
  }  
} # end of ploteHSphase




