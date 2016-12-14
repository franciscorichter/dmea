\name{llik_st}
\alias{llik_st}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{log-likelihood of a set of trees
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
llik_st(pars, setoftrees, impsam = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pars}{parameters the set of trees
%%     ~~Describe \code{pars} here~~
}
  \item{setoftrees}{ list with a set of trees
%%     ~~Describe \code{setoftrees} here~~
}
  \item{impsam}{ if impsam=TRUE, then the likelihood will be adjusted by importance scalings weights, included on the list of trees
%%     ~~Describe \code{impsam} here~~
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
function (pars, setoftrees, impsam = F)
{
    m = length(setoftrees)
    l = NULL
    for (i in 1:m) {
        s = setoftrees[[i]]
        if (impsam) {
            weight = s$weight
            if (weight == 0) {
                l[i] = 0
            }
            else {
                l[i] = llik(b = pars, n = s$n, E = s$E, t = s$wt) *
                  weight
            }
        }
        else {
            l[i] = llik(b = pars, n = s$n, E = s$E, t = s$wt)
        }
    }
    L = sum(l)
    return(L)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
