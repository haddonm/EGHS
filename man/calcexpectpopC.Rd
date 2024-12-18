% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tashs_funs.R
\name{calcexpectpopC}
\alias{calcexpectpopC}
\title{calcexpectpopC generates population catches from TAC or acatch}
\usage{
calcexpectpopC(
  hcrout,
  exb,
  sauCPUE,
  catsau,
  fleetacatch,
  hsargs = hsargs,
  glb = glb,
  sigmab = 0.1,
  year = year
)
}
\arguments{
\item{hcrout}{in Tasmania this contains the aspirational catches by SAU
along with many other details.}

\item{exb}{the exploitable biomass at the end of the previous year. In the
first year of the projections this would be the last year of the
conditioning.}

\item{sauCPUE}{the CPUE by sau, not used by TAS but needed by SA}

\item{catsau}{the actual catches by sau, not used by TAS but needed by SA}

\item{fleetacatch}{a function that calculates catch distribution across
the sau and populations. Currently not needed by Tas but needed by SA}

\item{hsargs}{the harvest strategy constants not used in TAS, default=hsargs}

\item{glb}{the globals object,  not used in TAS, default=glb}

\item{sigmab}{the Log-Normal standard deviation of implementation error. The
default value = 1e-08, which effectively means no errors.}

\item{year}{what year number are we dealing with?}
}
\value{
a vector of population catches for the year to be imposed after the
estimation of the mid-year exploitable biomass (midyexpB).
}
\description{
calcexpectpopC defines the fleet dynamics model used by the
operating model. The function converts the TAC into expected sau catches,
which are then translated into population level catches while
introducing management implementation error. Here this
error is implemented as Log-Normal errors on diver intuitions
concerning the relative abundance in each population. The error is
imposed separately on the populations in each SAU. The option of TAC for
a zone or aspirational catches for SAU is given to allow for different
processes between jurisdictions. If only a TAC is produced by an HCR then
calcexpectpopC needs to be able to generate catches by SAU prior to
subdividing those SAU catches among populations. Used directly by aMSE.
}
\examples{
print("wait on suitable internal data sets")
}
