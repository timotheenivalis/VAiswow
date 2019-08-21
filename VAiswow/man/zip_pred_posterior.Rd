\name{zip_pred_posterior}
\alias{zip_pred_posterior}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
zip_pred_posterior(m, nMC = 10, fitnessfile = NULL, fixedpred = c("inbreeding", "Qgg"), ranpred = c("id", "dam", "cohort"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{
%%     ~~Describe \code{m} here~~
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
function (m, nMC = 10, fitnessfile = NULL, fixedpred = c("inbreeding", 
    "Qgg"), ranpred = c("id", "dam", "cohort")) 
{
    length_posterior <- length(m$Sol[, 1])
    fit <- read.csv(fitnessfile)
    limdis <- max(fit$LBS, na.rm = TRUE) + 2
    predmat <- matrix(0, nrow = length_posterior, ncol = limdis + 
        1)
    for (i in 1:length_posterior) {
        pred <- zip_pred(m = m, s = i, nMC = nMC, fitnessfile = fitnessfile, 
            fixedpred = fixedpred, ranpred = ranpred)
        pred[pred > limdis] <- limdis
        tabp <- table(pred)
        tabp <- tabp/sum(tabp)
        tabpnames <- as.numeric(names(tabp))
        for (j in 1:length(tabpnames)) {
            predmat[i, tabpnames[j] + 1] <- tabp[j]
        }
    }
    return(predmat)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
