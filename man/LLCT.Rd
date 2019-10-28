\name{LLCT}
\alias{LLCT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Longitudinal Linear Combination Test for Gene Set Analysis
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
  LLCT is a two-step self-contained gene-set analysis method which is developed to handle multiple longitudinal outcomes. Analysis of within-subject variation in the first step is followed by examining the between-subject variation utilizing Linear Combination Test (LCT) in the second step. This method is also applicable in analysis of time-course microarray data.
}
\usage{
LLCT(EXPR, GS, LongData, ID = "ID", time = "time", covariate = NULL, phenotype = "phenotype", familybased = FALSE, pedigree = NULL, FIX.formula = "phenotype~time+covariate", RANDOM.formula = "~1|ID", nbPermutations = 1000, family = "gaussian(link=identity)")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{EXPR}{
%%     ~~Describe \code{EXPR} here~~
Gene Expressions --- Row names:Subjects (single row per subject); Col names: Genes names. Note: The rownames of GS must be matched with colnames of EXPR and the unique values of the rownames of LongData must be matched with the rownames of EXPR.   
}
  \item{GS}{
%%     ~~Describe \code{GS} here~~
    Gene Set Matrix --- Row names: list of genes; Col names: list of selected gene sets; Cells: 1 (if gene of the row belongs to the gene set of the column) 0 (otherwise). Note: The rownames of GS must be matched with colnames of EXPR
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
data(data_LLCT)

#LLCT for unrelated subjects with two time variables (time and time2) defining phenotypic temppral patterns 
LLCT(EXPR=EXPR, GS=GS, LongData=LongData, ID="ID",time=c("time","time2"), covariate=NULL,phenotype="phenotype",familybased=FALSE,pedigree=NULL,FIX.formula="phenotype~time+time2", RANDOM.formula=NULL,nbPermutations=1000)

#LLCT for unrelated subjects with a linear time trend 
#LLCT(EXPR=EXPR, GS=GS, LongData=LongData, ID="ID",time=c("time"), covariate=NULL,phenotype="phenotype",familybased=FALSE,pedigree=NULL,FIX.formula="phenotype~time", RANDOM.formula=NULL,nbPermutations=1000)

#LLCT for unrelated subjects with binary phenotype 
#LLCT(EXPR=EXPR, GS=GS, LongData=LongData, ID="ID",time=c("time"), covariate=NULL,phenotype="binary.phenotype",familybased=FALSE,pedigree=NULL,FIX.formula="binary.phenotype~time", RANDOM.formula=NULL,nbPermutations=1000,family="binomial(link=logit)")

#LLCT for unrelated subjects with adjustment for time dependent covariate
#LLCT(EXPR=EXPR, GS=GS, LongData=LongData, ID="ID",time=c("time"), covariate="covariate",phenotype="phenotype",familybased=TRUE,pedigree="pedigree",FIX.formula="phenotype~time+covariate", RANDOM.formula="~1|ID",nbPermutations=1000)

#LLCT for related subjects with adjustment for time dependent covariate 
#LLCT(EXPR=EXPR, GS=GS, LongData=LongData, ID="ID",time=c("time"), covariate="covariate",phenotype="phenotype",familybased=TRUE,pedigree="pedigree",FIX.formula="phenotype~time+covariate", RANDOM.formula="~1|ID",nbPermutations=1000)

#LLCT for related subjects with adjustment for more complex time trend 
#LLCT(EXPR=EXPR, GS=GS, LongData=LongData, ID="ID",time=c("time","time2"), covariate="covariate",phenotype="phenotype",familybased=TRUE,pedigree="pedigree",FIX.formula="phenotype~time+time2", RANDOM.formula="~1|ID",nbPermutations=1000)

#LLCT for related subjects with binary phenotype 
#LLCT(EXPR=EXPR, GS=GS, LongData=LongData, ID="ID",time=c("time"), covariate="covariate",phenotype="binary.phenotype",familybased=TRUE,pedigree="pedigree",FIX.formula="binary.phenotype~time", RANDOM.formula="~1|ID",nbPermutations=1000,family="binomial(link=logit)")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
