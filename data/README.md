After loadind data_LLCT.rda:

#LLCT for unrelated subjects with two time variables (time and time2) defining phenotypic temppral patterns
LLCT(EXPR=EXPR, GS=GS, LongData=LongData, ID="ID",time=c("time","time2"), covariate=NULL,phenotype="phenotype",familybased=FALSE,pedigree=NULL,FIX.formula="phenotype~time+time2", RANDOM.formula=NULL,nbPermutations=1000)

#LLCT for unrelated subjects with a linear time trend
LLCT(EXPR=EXPR, GS=GS, LongData=LongData, ID="ID",time=c("time"), covariate=NULL,phenotype="phenotype",familybased=FALSE,pedigree=NULL,FIX.formula="phenotype~time", RANDOM.formula=NULL,nbPermutations=1000)

#LLCT for unrelated subjects with binary phenotype
LLCT(EXPR=EXPR, GS=GS, LongData=LongData, ID="ID",time=c("time"), covariate=NULL,phenotype="binary.phenotype",familybased=FALSE,pedigree=NULL,FIX.formula="binary.phenotype~time", RANDOM.formula=NULL,nbPermutations=1000,family="binomial(link=logit)")

#LLCT for unrelated subjects with adjustment for time dependent covariate
LLCT(EXPR=EXPR, GS=GS, LongData=LongData, ID="ID",time=c("time"), covariate="covariate",phenotype="phenotype",familybased=TRUE,pedigree="pedigree",FIX.formula="phenotype~time+covariate", RANDOM.formula="~1|ID",nbPermutations=1000)

#LLCT for related subjects with adjustment for time dependent covariate
LLCT(EXPR=EXPR, GS=GS, LongData=LongData, ID="ID",time=c("time"), covariate="covariate",phenotype="phenotype",familybased=TRUE,pedigree="pedigree",FIX.formula="phenotype~time+covariate", RANDOM.formula="~1|ID",nbPermutations=1000)

#LLCT for related subjects with adjustment for more complex time trend
LLCT(EXPR=EXPR, GS=GS, LongData=LongData, ID="ID",time=c("time","time2"), covariate="covariate",phenotype="phenotype",familybased=TRUE,pedigree="pedigree",FIX.formula="phenotype~time+time2", RANDOM.formula="~1|ID",nbPermutations=1000)

#LLCT for related subjects with binary phenotype
LLCT(EXPR=EXPR, GS=GS, LongData=LongData, ID="ID",time=c("time"), covariate="covariate",phenotype="binary.phenotype",familybased=TRUE,pedigree="pedigree",FIX.formula="binary.phenotype~time", RANDOM.formula="~1|ID",nbPermutations=1000,family="binomial(link=logit)")
