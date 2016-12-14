\name{phylo2p}
\alias{phylo2p}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Convert phylo object on a list with 3 vectors: waiting times, number of species at each node and a binary vector with  speciation/extinction identification for each node.
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
phylo2p(tree)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{Object of class phylo
%%     ~~Describe \code{tree} here~~
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
function (tree)
{
    ltt = ltt.plot.coords(tree)
    t = diff(ltt[, 1])
    ltt = ltt[-1, ]
    n = ltt[, 2]
    E = diff(n)
    n = n[-length(n)]
    E[E == -1] = 0
    return(list(t = t, E = E, n = n))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
