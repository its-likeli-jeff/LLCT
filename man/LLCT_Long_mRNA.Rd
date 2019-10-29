\name{LLCT_Long_mRNA}
\alias{LLCT_Long_mRNA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ 
  Longitudinal Linear Combination Test for Gene Set Analysis
}
\description{
  LLCT is a two-step self-contained gene-set analysis method which is developed to handle multiple longitudinal outcomes. Analysis of within-subject variation in the first step is followed by examining the between-subject variation utilizing Linear Combination Test (LCT) in the second step. This method is also applicable in analysis of time-course microarray data.
}
\usage{
LLCT_Long_mRNA(LongEXPR=LongEXPR, GS=GS, PhenData=PhenData, ID="ID",time=c("time","time2"), covariate=NULL, Genes=colnames(LongEXPR)[!(colnames(LongEXPR) %in% c("ID","time","time2","covariate"))],phenotype="phenotype",FIX.formula="~time+time2",nbPermutations=1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{LongEXPR}{
Gene Expressions --- Row names:Subjects (multiple rows per subject); Col names: columns include Subjects' ID, time/visit variables, time dependent covariates and Genes expressions . The rownames of GS must be matched with  names of the genes included in LongEXPR. The ID values of PhenData must be matched with unique rownames of LongEXPR.
}
  \item{GS}{
    Gene Set Matrix --- Row names: list of genes; Col names: list of selected gene sets; Cells: 1 (if gene of the row belongs to the gene set of the column) 0 (otherwise). Note: The rownames of GS must be matched with colnames of EXPR
}
  \item{PhenData}{
Phenotype(s)/Predictor(s) Data. Columns should include Subjects' ID, single or multiple phenotypes. Rows: Subjects (one row per subject). The ID values of PhenData must be matched with the unique rownames of LongEXPR 
}
  \item{ID}{
Name of ID variable in PhenData and LongEXPR Defaults to "ID".
}
  \item{time}{
Vector of the names of time variables in LongEXPR
}
  \item{covariate}{
name of covariate(s) in LongEXPR
}
\item{Genes}{
List of gene names within LongEXPR defaults to colnames(LongEXPR)[-1:-2].
}
  \item{phenotype}{
 name of phenotype in PhenData
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
data(data_for_LLCT_Long_mRNA)
#LLCT for time-course microarray datawith two time variables (time and time2) defining phenotypic temppral patterns
mrnarun <- LLCT_Long_mRNA(LongEXPR=LongEXPR, GS=GS, PhenData=PhenData, ID="ID",time=c("time","time2"), covariate=NULL, Genes=colnames(LongEXPR)[!(colnames(LongEXPR) %in% c("ID","time","time2","covariate"))],phenotype="phenotype",FIX.formula="~time+time2",nbPermutations=1000)
mrnarun$LLCT_Results
  
#'#LLCT for time-course microarray data with a linear time trend
#'LLCT_Long_mRNA(LongEXPR=LongEXPR, GS=GS, PhenData=PhenData, ID="ID",time=c("time"), covariate=NULL, Genes=colnames(LongEXPR)[!(colnames(LongEXPR) %in% c("ID","time","time2","covariate"))],phenotype="phenotype",FIX.formula="~time+time2",nbPermutations=1000)


#'#LLCT for time-course microarray data with adjustment for time dependent covariate
#'LLCT_Long_mRNA(LongEXPR=LongEXPR, GS=GS, PhenData=PhenData, ID="ID",time=c("time"), covariate="covariate", Genes=colnames(LongEXPR)[!(colnames(LongEXPR) %in% c("ID","time","time2","covariate"))],phenotype="phenotype",FIX.formula="~time+time2",nbPermutations=1000)

#'#LLCT for time-course microarray data with two phenotypes
#'LLCT_Long_mRNA(LongEXPR=LongEXPR, GS=GS, PhenData=PhenData, ID="ID",time=c("time"), covariate=NULL, Genes=colnames(LongEXPR)[!(colnames(LongEXPR) %in% c("ID","time","time2","covariate"))],phenotype=c("phenotype","phenotype2"),FIX.formula="~time+time2",nbPermutations=1000)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
