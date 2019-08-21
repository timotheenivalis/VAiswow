\name{Va_w_zip_posterior}
\alias{Va_w_zip_posterior}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Va_w_zip_posterior(m, nMC = 50000, h = 0.02, fitnessfile, fixedpred = c("inbreeding", "Qgg"), ranpred = c("id", "dam", "cohort"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{
%%     ~~Describe \code{m} here~~
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
function (m, nMC = 50000, h = 0.02, fitnessfile, fixedpred = c("inbreeding", 
    "Qgg"), ranpred = c("id", "dam", "cohort")) 
{
    length_posterior <- length(m$Sol[, 1])
    ranpredunits <- c(ranpred, "units")
    btzip <- as.mcmc(matrix(data = NA, nrow = length_posterior, 
        ncol = length(ranpredunits) + 1, dimnames = list(1:length_posterior, 
            c(ranpredunits, "residual"))))
    attr(btzip, "mcpar") <- attr(m$VCV, "mcpar")
    for (va in 1:length(ranpredunits)) {
        btva <- vector(length = length_posterior)
        for (i in 1:length_posterior) {
            btva[i] <- Va_w_zip(m = m, s = i, fitnessfile = fitnessfile, 
                fixedpred = fixedpred, nMC = nMC, h = h, ranpred = ranpred, 
                focalvar = ranpredunits[va])
        }
        btzip[, ranpredunits[va]] <- btva
    }
    fit <- read.csv(fitnessfile)
    btzip[, "residual"] <- var(fit$LBS) - rowSums(btzip[, -ncol(btzip)])
    m$btzip <- btzip
    return(m)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
