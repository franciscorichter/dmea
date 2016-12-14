\name{EM_phylo}
\alias{EM_phylo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A montecarlo EM algorithm implementation to diversification rates inference
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
EM_phylo(wt, init_par, n_trees = 100, n_it = 30, printpar = TRUE, mu = NULL, impsam = FALSE, dummy = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{wt}{ waiting times of the (extant only) phylogenetic tree
%%     ~~Describe \code{wt} here~~
}
  \item{init_par}{ Initial parameters for the EM algorithm
%%     ~~Describe \code{init_par} here~~
}
  \item{n_trees}{ Number of trees to be samplead for the montecarlo algorithm at each EM iteration
%%     ~~Describe \code{n_trees} here~~
}
  \item{n_it}{ Number of iterations for the EM algorithm (this should be replaced by a tolerance value and a stoping criteria)
%%     ~~Describe \code{n_it} here~~
}
  \item{printpar}{if printpar=TRUE then the EM estimations will be printed at each iteration
%%     ~~Describe \code{printpar} here~~
}
  \item{mu}{ You can use this in case you do not want to estimate mu
%%     ~~Describe \code{mu} here~~
}
  \item{impsam}{ if impsam=TRUE then an importance sampling approach will be used on the montecarlo sampling
%%     ~~Describe \code{impsam} here~~
}
  \item{dummy}{
%%     ~~Describe \code{dummy} here~~
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
function (wt, init_par, n_trees = 100, n_it = 30, printpar = TRUE,
    mu = NULL, impsam = FALSE, dummy = 0)
{
    pars = init_par
    for (j in 1:n_it) {
        if (printpar)
            print(pars)
        trees <- sim_srt(wt = wt, pars = pars, parallel = F,
            n_trees = n_trees)
        if (length(init_par) == 3) {
            pars = subplex(par = c(8, 0.175, 0.9), fn = llik_st,
                setoftrees = trees, impsam = impsam)$par
        }
        else {
            pars = subplex(par = c(8, 0.175), mu2par = mu, fn = llik_st,
                setoftrees = trees, impsam = impsam)$par
        }
        pars[3] = pars[3] + dummy
    }
    return(pars)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
