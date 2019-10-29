#' LLCT Function for Time-Course Microarray Data
#'
#' This function allows you to conduct a time-course gene set enrichment anlysis 
#' LLCT in this application tests if the gene sets are differentially expressed over time in association with an underlying phenotype(s) or predictor(s)
#' @param GS Gene Set Matrix
#' Row names: list of genes
#' Col names: list of selected gene sets
#' Cells: 1 (if gene of the row belongs to the gene set of the column) 0 (otherwise)
#' The rownames of GS must be matched with colnames of EXPR
#' @param PhenData Phenotype(s)/Predictor(s) Data
#' Columns should include Subjects' ID, single or multiple phenotypes 
#' Rows: Subjects (one row per subject)
#' The ID values of PhenData must be matched with the unique rownames of LongEXPR 
#' @param LongEXPR Gene Expressions
#' Row names:Subjects (multiple rows per subject)
#' Col names: columns include Subjects' ID, time/visit variables, time dependent covariates amd Genes expressions
#' The rownames of GS must be matched with  names of the genes included in LongEXPR 
#' The ID values of PhenData must be matched with unique rownames of LongEXPR 
#' @param Genes List of gene names within LongEXPR defaults to colnames(LongEXPR)[-1:-2].
#' @param ID Name of ID variable in PhenData and LongEXPR Defaults to "ID".
#' @param time vector of the names of time variables in LongEXPR Defaults to "time".
#' @param covariate name of covariate(s) in LongEXPR Defaults to NULL.
#' @param phenotype Defaults to  names of phenotypes in PhenData Default to "phenotype"
#' @param FIX.formula Formula of Phenotype~time model Defaults to "~time+covariate"
#' @param nbPermutations Number of Permutations Defaults to 1000
#' @return "LLCT_Results" geneset-specific p-values and q-values
#' @return "Step1_Coefs" coefficients calculated in the first step using FIX_formula 
#' @keywords LLCT
#' @export 
#' @examples
#'#LLCT for time-course microarray datawith two time variables (time and time2) defining phenotypic temppral patterns
#'LLCT_Long_mRNA(LongEXPR=LongEXPR, GS=GS, PhenData=PhenData, ID="ID",time=c("time","time2"), covariate=NULL, Genes=colnames(LongEXPR)[!(colnames(LongEXPR) %in% c("ID","time","time2","covariate"))],phenotype="phenotype",FIX.formula="~time+time2",nbPermutations=1000)

#'#LLCT for time-course microarray data with a linear time trend
#'LLCT_Long_mRNA(LongEXPR=LongEXPR, GS=GS, PhenData=PhenData, ID="ID",time=c("time"), covariate=NULL, Genes=colnames(LongEXPR)[!(colnames(LongEXPR) %in% c("ID","time","time2","covariate"))],phenotype="phenotype",FIX.formula="~time+time2",nbPermutations=1000)


#'#LLCT for time-course microarray data with adjustment for time dependent covariate
#'LLCT_Long_mRNA(LongEXPR=LongEXPR, GS=GS, PhenData=PhenData, ID="ID",time=c("time"), covariate="covariate", Genes=colnames(LongEXPR)[!(colnames(LongEXPR) %in% c("ID","time","time2","covariate"))],phenotype="phenotype",FIX.formula="~time+time2",nbPermutations=1000)

#'#LLCT for time-course microarray data with two phenotypes
#'LLCT_Long_mRNA(LongEXPR=LongEXPR, GS=GS, PhenData=PhenData, ID="ID",time=c("time"), covariate=NULL, Genes=colnames(LongEXPR)[!(colnames(LongEXPR) %in% c("ID","time","time2","covariate"))],phenotype=c("phenotype","phenotype2"),FIX.formula="~time+time2",nbPermutations=1000)



LLCT_Long_mRNA=function(LongEXPR, GS, PhenData, ID="ID", time="time", covariate=NULL, phenotype="phenotype", Genes=colnames(LongEXPR)[-1:-2],  FIX.formula="~time+covariate", nbPermutations=1000){




	nb.Samples<-length(PhenData$ID) #Number of Subjects
	nb.Genes<-length(Genes) #Number of Genes
	model.coefs<-NULL

	#Ordering all data
	GS<-GS[,order(colnames(GS))]
	GS<-GS[order(rownames(GS)),]

	

		#Step 1 for unrelated subjects:
		betas3=NULL
		for (k in 1:nb.Genes){
		model.coefs<-lapply(split(LongEXPR,LongEXPR$ID),function(z){as.data.frame(summary(lm(formula(paste("z[,colnames(z)==Genes[k]]",FIX.formula)),data=z))$coefficients)})
										#Calculation of Subject specific trends
		betas<-lapply(time,function(w){lapply(model.coefs,function(z){(z)[which(rownames(z) %in% w ),"Estimate"]})})
		betas2<-matrix(unlist(lapply(betas,function(z)lapply(z,function(w){if(length(w)==0){w=NA}else{w}}))),ncol=length(time))
		if(sum(is.na(betas2))>0){print("The formula did not fit all subjects")}
		betas3=cbind(betas2,betas3)
		}

		
		rownames(betas3)=names(model.coefs)
		#colnames(betas3)=Genes
		betas3=betas3[order(rownames(betas3)),]
		#betas3=betas3[,order(colnames(betas3))]



	#Step 2: Analyses of Between Subject Variations
	GS=GS[order(rownames(GS)),]
	PhenData=PhenData[order(PhenData$ID),]
	GS2= as.data.frame((GS))
	GS.sizes <- sapply(GS2,sum) # size of each gene set
	GS.data<-NULL	
	GS.data <- lapply(data.frame(GS2), function(z) (betas3)[,z==1]); # creat data of each GS (rows=gene trends, columns=samples/families)
	#GS.data <- lapply(GS.data,function(z) t(z)); # transfer genes in each GS (columns=gene trends, rows=samples/families)
	# (2) Eigen-decompsition of shrinkage pooled covariance matrix for each GS
	Cov.Pooled<-lapply(GS.data, function(z) cov.shrink(z,verbose=FALSE, lambda.var=0));   # pooled covariance of genes in each GS

	nb.GeneSets<-dim(GS)[2]
                          
   	for (i in 1:nb.GeneSets){
         EIGEN.decom<-eigen(Cov.Pooled[[i]]);
                                     # eigen decomposition of pooled covariance for each GS
         D<-EIGEN.decom$values;      
         U<-EIGEN.decom$vectors;
		if(sum(D<0)>0){
			class(Cov.Pooled[[i]])="matrix"
			EIGEN.decom<-eigen(nearPD((Cov.Pooled[[i]]))[[1]]);
			D<-EIGEN.decom$values;      
      		U<-EIGEN.decom$vectors;

		}
         GS.data[[i]]<-t(GS.data[[i]]%*%U)/sqrt(D)
                                     # adjust data of each GS (rows=genes, columns=samples)
      }

	PHEN=as.matrix(PhenData[,phenotype])
      Cov.Pooled<-cov.shrink(PHEN,verbose=FALSE, lambda.var=0);
      EIGEN.decom<-eigen(Cov.Pooled);  # eigen decomposition of pooled covariance for beta

	D<-EIGEN.decom$values;         
      U<-EIGEN.decom$vectors;

      cl<-t(PHEN%*%U)/sqrt(D)          # adjust beta (rows=respons, columns=samples)
	
	
	# (3) T-like stats obtained on 'true' data
      sam.sumsquareT.obs  <- sapply(GS.data, function(z) T2.like.SAMGS(z,cl))
                                      # the T-like statistics obtained on 'true' data


     # (4) stats obtained on 'permuted' data
     sam.sumsquareT.permut <- matrix(NA,nbPermutations,nb.GeneSets)
     for(i in 1:nbPermutations) {

		
      		ind <- sample(nb.Samples)
		
       sam.sumsquareT.permut[i,] <- sapply(GS.data, function(z) T2.like.SAMGS(z[,ind],cl))
                                      # SAMGS statitic for each gene set  - for current permutation
      }

	# (5) p-value and q-value
	GeneSets.pval <- apply(t(sam.sumsquareT.permut) >= sam.sumsquareT.obs,1,sum)/nbPermutations

	GeneSets.qval <-0; 
	GeneSets.qval <-p.adjust(GeneSets.pval,method="fdr")

	res <- as.data.frame(cbind("GS size"              = GS.sizes,
                                "GS p-value" 	       = GeneSets.pval
                               ,"GS q-value"           = GeneSets.qval))




	return(list("LLCT_Results"=res, "Step1_Coefs"=betas3))

}



