\name{sim_srt}
\alias{sim_srt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulate a set of complete trees given the observed tree
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sim_srt(wt, pars, parallel = F, n_trees)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{wt}{waiting times of extant species
%%     ~~Describe \code{wt} here~~
}
  \item{pars}{parameters for the reconstruction
%%     ~~Describe \code{pars} here~~
}
  \item{parallel}{if parallel=TRUE then the algorithm will be running in parallel using n-1 of the cores of your computer
%%     ~~Describe \code{parallel} here~~
}
  \item{n_trees}{number of trees to simulate
%%     ~~Describe \code{n_trees} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{list with n_trees simulated complete trees
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
function (wt, pars, parallel = F, n_trees)
{
    if (parallel) {
        trees <- foreach(i = 1:n_trees, combine = list) \%dopar\%
            dmea::rec_tree(w = tree$t, pars = pars)
    }
    else {
        trees = vector("list", length = n_trees)
        for (i in 1:n_trees) {
            rec = dmea::rec_tree(wt = wt, pars = pars)
            trees[[i]] = rec
        }
    }
    return(trees)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
