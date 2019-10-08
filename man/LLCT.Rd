\name{LLCT}
\alias{LLCT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
LLCT(EXPR, GS, LongData, ID = "ID", time = "time", covariate = NULL, phenotype = "phenotype", familybased = FALSE, pedigree = NULL, FIX.formula = "phenotype~time+covariate", RANDOM.formula = "~1|ID", nbPermutations = 1000, family = "gaussian(link=identity)")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{EXPR}{
%%     ~~Describe \code{EXPR} here~~
}
  \item{GS}{
%%     ~~Describe \code{GS} here~~
}
  \item{LongData}{
%%     ~~Describe \code{LongData} here~~
}
  \item{ID}{
%%     ~~Describe \code{ID} here~~
}
  \item{time}{
%%     ~~Describe \code{time} here~~
}
  \item{covariate}{
%%     ~~Describe \code{covariate} here~~
}
  \item{phenotype}{
%%     ~~Describe \code{phenotype} here~~
}
  \item{familybased}{
%%     ~~Describe \code{familybased} here~~
}
  \item{pedigree}{
%%     ~~Describe \code{pedigree} here~~
}
  \item{FIX.formula}{
%%     ~~Describe \code{FIX.formula} here~~
}
  \item{RANDOM.formula}{
%%     ~~Describe \code{RANDOM.formula} here~~
}
  \item{nbPermutations}{
%%     ~~Describe \code{nbPermutations} here~~
}
  \item{family}{
%%     ~~Describe \code{family} here~~
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

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
