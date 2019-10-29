\name{LLCT}
\alias{LLCT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ 
  Longitudinal Linear Combination Test for Gene Set Analysis
}
\description{
  LLCT is a two-step self-contained gene-set analysis method which is developed to handle multiple longitudinal outcomes. Analysis of within-subject variation in the first step is followed by examining the between-subject variation utilizing Linear Combination Test (LCT) in the second step. This method is also applicable in analysis of time-course microarray data.
}
\usage{
LLCT(EXPR, GS, LongData, ID = "ID", time = "time", covariate = NULL, phenotype = "phenotype", familybased = FALSE, pedigree = NULL, FIX.formula = "phenotype~time+covariate", RANDOM.formula = "~1|ID", nbPermutations = 1000, family = "gaussian(link=identity)")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{EXPR}{
Gene Expressions --- Row names:Subjects (single row per subject); Col names: Genes names. Note: The rownames of GS must be matched with colnames of EXPR and the unique values of the rownames of LongData must be matched with the rownames of EXPR.   
}
  \item{GS}{
    Gene Set Matrix --- Row names: list of genes; Col names: list of selected gene sets; Cells: 1 (if gene of the row belongs to the gene set of the column) 0 (otherwise). Note: The rownames of GS must be matched with colnames of EXPR
}
  \item{LongData}{
Phenotypes and Covariates Data in Longitudinal format
}
  \item{ID}{
Name of ID variable in LongData
}
  \item{time}{
Vector of the names of time variables in LongData
}
  \item{covariate}{
name of covariate(s) in LongData
}
  \item{phenotype}{
 name of phenotype in LongData
}
  \item{familybased}{
TRUE/FALSE indicating if the subjects are sharing families/pedigress Defaults to FALSE
}
  \item{pedigree}{
name of pedigree variable in LongData
}
  \item{FIX.formula}{
Formula of Phenotype~time model Defaults to "phenotype~time+covariate"
}
  \item{RANDOM.formula}{
Formula of random effects as lme function requires - required in case of related subjects Defaults to NULL
}
  \item{nbPermutations}{
Number of Permutations Defaults to 1000
}
  \item{family}{
family of link function appropriate for the type of phenotype
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{LLCT_Results}{geneset-specific p-values and q-values}
\item{Step1_Coefs}{coefficients calculated in the first step using FIX_formula}
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

#LLCT for unrelated subjects with two more complex time trend
LLCT(EXPR=EXPR, GS=GS, LongData=LongData, ID="ID",time=c("time","time2"), covariate=NULL,phenotype="phenotype",familybased=FALSE,pedigree=NULL,FIX.formula="~time+time2", RANDOM.formula=NULL,nbPermutations=1000)

#LLCT for unrelated subjects with a linear time trend
#LLCT(EXPR=EXPR, GS=GS, LongData=LongData, ID="ID",time=c("time"), covariate=NULL,phenotype="phenotype",familybased=FALSE,pedigree=NULL,FIX.formula="~time", RANDOM.formula=NULL,nbPermutations=1000)

#LLCT for unrelated subjects with adjustment for time dependent covariate
#LLCT(EXPR=EXPR, GS=GS, LongData=LongData, ID="ID",time=c("time"), covariate="covariate",phenotype="phenotype",familybased=FALSE,pedigree=NULL,FIX.formula="~time+covariate", RANDOM.formula=NULL,nbPermutations=1000)

#LLCT for unrelated subjects with binary phenotype
#LLCT(EXPR=EXPR, GS=GS, LongData=LongData, ID="ID",time=c("time"), covariate=NULL,phenotype="binary.phenotype",familybased=FALSE,pedigree=NULL,FIX.formula="~time", RANDOM.formula=NULL,nbPermutations=1000,family="binomial(link=logit)")

#LLCT for related subjects with adjustment for time dependent covariate
#LLCT(EXPR=EXPR, GS=GS, LongData=LongData, ID="ID",time=c("time"), covariate="covariate",phenotype="phenotype",familybased=FALSE,pedigree="pedigree",FIX.formula="~time+covariate", RANDOM.formula="~1|ID",nbPermutations=1000)

#LLCT for related subjects with adjustment for more complex time trend
#LLCT(EXPR=EXPR, GS=GS, LongData=LongData, ID="ID",time=c("time","time2"), covariate="covariate",phenotype="phenotype",familybased=FALSE,pedigree="pedigree",FIX.formula="~time+time2", RANDOM.formula="~1|ID",nbPermutations=1000)

#LLCT for unrelated subjects with two phenotypes
#LLCT(EXPR=EXPR, GS=GS, LongData=LongData, ID="ID",time=c("time"), covariate=NULL,phenotype=c("phenotype","phenotype2"),familybased=FALSE,pedigree=NULL,FIX.formula="~time", RANDOM.formula=NULL,nbPermutations=1000)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
