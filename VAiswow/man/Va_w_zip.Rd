\name{Va_w_zip}
\alias{Va_w_zip}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Va_w_zip(m, s = 1, nMC = 50000, h = 0.02, fitnessfile = NULL, fixedpred = c("inbreeding", "Qgg"), ranpred = c("id", "dam", "cohort"), focalvar = "id")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{
%%     ~~Describe \code{m} here~~
}
  \item{s}{
%%     ~~Describe \code{s} here~~
}
  \item{nMC}{
%%     ~~Describe \code{nMC} here~~
}
  \item{h}{
%%     ~~Describe \code{h} here~~
}
  \item{fitnessfile}{
%%     ~~Describe \code{fitnessfile} here~~
}
  \item{fixedpred}{
%%     ~~Describe \code{fixedpred} here~~
}
  \item{ranpred}{
%%     ~~Describe \code{ranpred} here~~
}
  \item{focalvar}{
%%     ~~Describe \code{focalvar} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (m, s = 1, nMC = 50000, h = 0.02, fitnessfile = NULL, 
    fixedpred = c("inbreeding", "Qgg"), ranpred = c("id", "dam", 
        "cohort"), focalvar = "id") 
{
    muP <- m$Sol[s, "traitLBS"]
    muZI <- m$Sol[s, "traitzi_LBS"]
    if (!is.null(fixedpred)) {
        if (is.null(fitnessfile)) {
            stop("No fitnessfile provided but need to integrate over fixed effects")
        }
        else {
            fit <- read.csv(fitnessfile)
        }
        for (i in 1:length(fixedpred)) {
            muP <- muP + fit[, fixedpred[i]] * m$Sol[s, paste0("traitLBS:", 
                fixedpred[i])]
            muZI <- muZI + fit[, fixedpred[i]] * m$Sol[s, paste0("traitzi_LBS:", 
                fixedpred[i])]
        }
    }
    mu <- cbind(muP, muZI)
    PVP <- m$VCV[s, "traitLBS.units"] + var(mu[, 1])
    for (i in 1:length(ranpred)) {
        PVP <- PVP + m$VCV[s, paste0("traitLBS:traitLBS.", ranpred[i])]
    }
    PVZI <- m$VCV[s, "traitzi_LBS.units"] + var(mu[, 2])
    for (i in 1:length(ranpred)) {
        PVZI <- PVZI + m$VCV[s, paste0("traitzi_LBS:traitzi_LBS.", 
            ranpred[i])]
    }
    PcovZIP <- cov(mu[, 1], mu[, 2])
    for (i in 1:length(ranpred)) {
        PcovZIP <- PcovZIP + m$VCV[s, paste0("traitzi_LBS:traitLBS.", 
            ranpred[i])]
    }
    sigma <- matrix(c(PVP, PcovZIP, PcovZIP, PVZI), 2, 2)
    MC_latent_smaples <- rmvnorm(nMC, colMeans(mu), sigma)
    muW <- barW(MC_latent_smaples)
    muW1 <- barW(MC_latent_smaples + cbind(rep(h, nMC), 0))
    muW2 <- barW(MC_latent_smaples + cbind(0, rep(h, nMC)))
    d1 <- (muW1/muW - 1)/h
    d2 <- (muW2/muW - 1)/h
    if (focalvar == "units") {
        sigma_a <- matrix(data = c(m$VCV[s, paste0("traitLBS.", 
            focalvar)], 0, 0, m$VCV[s, paste0("traitzi_LBS.", 
            focalvar)]), 2, 2)
    }
    else {
        sigma_a <- matrix(m$VCV[s, c(paste0(c("traitLBS:traitLBS.", 
            "traitzi_LBS:traitLBS.", "traitLBS:traitzi_LBS.", 
            "traitzi_LBS:traitzi_LBS."), focalvar))], 2, 2)
    }
    d <- matrix(c(d1, d2), 2, 1)
    return(as.numeric(t(d) \%*\% sigma_a \%*\% d))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
