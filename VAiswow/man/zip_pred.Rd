\name{zip_pred}
\alias{zip_pred}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
zip_pred(m, s = 1, nMC = 50, fitnessfile = NULL, fixedpred = c("inbreeding", "Qgg"), ranpred = c("id", "dam", "cohort"))
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
  \item{fitnessfile}{
%%     ~~Describe \code{fitnessfile} here~~
}
  \item{fixedpred}{
%%     ~~Describe \code{fixedpred} here~~
}
  \item{ranpred}{
%%     ~~Describe \code{ranpred} here~~
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
function (m, s = 1, nMC = 50, fitnessfile = NULL, fixedpred = c("inbreeding", 
    "Qgg"), ranpred = c("id", "dam", "cohort")) 
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
    pred <- vector(length = nrow(mu) * nMC)
    for (i in 1:nrow(mu)) {
        MC_latent_smaples <- rmvnorm(nMC, mu[i, ], sigma)
        pred0 <- sapply(MC_latent_smaples[, 2], function(x) {
            rbinom(n = 1, size = 1, prob = 1/(1 + exp(x)))
        })
        predcond <- sapply(MC_latent_smaples[, 1], function(x) {
            rpois(n = 1, lambda = exp(x))
        })
        pred[(1 + (i - 1) * (nMC)):(nMC * i)] <- pred0 * predcond
    }
    return(pred)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
