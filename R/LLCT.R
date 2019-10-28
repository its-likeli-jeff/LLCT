LLCT=function(EXPR, GS, LongData, ID="ID", time="time", covariate=NULL, phenotype="phenotype", familybased=FALSE, pedigree=NULL, FIX.formula="~time+covariate", RANDOM.formula="~1|ID", nbPermutations=1000, family="gaussian(link=identity)"){

	nb.Samples<-length(unique(LongData$ID)) #Number of Subjects
	model.coefs<-NULL

	#Ordering all data
	EXPR<- EXPR[order(rownames(EXPR)),]
	EXPR<- EXPR[,order(colnames(EXPR))]
	GS<-GS[,order(colnames(GS))]


	if (familybased==FALSE){ #Analysis of Unrelated Subjects

		#Step 1 for unrelated subjects:
		betas3=NULL
		for (k in 1:length(phenotype)){
		model.coefs<-lapply(split(LongData,LongData$ID),function(z){as.data.frame(summary(glm(formula(paste(phenotype[k],FIX.formula)),data=z,family=eval(parse(text =family))))$coefficients)})
										#Calculation of Subject specific trends
		betas<-lapply(time,function(w){lapply(model.coefs,function(z){(z)[which(rownames(z) %in% w ),"Estimate"]})})
		betas2<-matrix(unlist(lapply(betas,function(z)lapply(z,function(w){if(length(w)==0){w=NA}else{w}}))),ncol=length(time))
		if(sum(is.na(betas2))>0){print("The formula did not fit all subjects")}
		betas3=cbind(betas2,betas3)
		}

		EXPR2<-t(EXPR)
		rownames(betas3)=names(model.coefs)
	}else{ #Analysis of Related Subjects
		#Step 1 for related subjects - Fitting Phenotype:
		nb.Families<-length(unique(LongData$pedigree))

		
		if (family=="gaussian(link=identity)"){ #Analysis of Continuous Phenotype
		betas3=NULL
			for(p in 1:length(phenotype)){
				model.coefs<-lapply(split(LongData,LongData$pedigree),function(z){as.data.frame(summary(lme(formula(paste(phenotype[p],FIX.formula)),random=formula(RANDOM.formula),data=z))$coefficients$fixed)})
										#Calculation of Family specific trends	
				betas<-lapply(time,function(w){lapply(model.coefs,function(z){(z)[which(rownames(z) %in% w ),]})})
				betas2<-matrix(unlist(lapply(betas,function(z)lapply(z,function(w){if(length(w)==0){w=NA}else{w}}))),ncol=length(time))
				betas3=cbind(betas2,betas3)
			}
		}else{	#Analysis of binary/categorical Phenotypes
		betas3=NULL					
			for(p in 1:length(phenotype)){
				model.coefs<-lapply(split(LongData,LongData$pedigree),function(z){as.data.frame(summary(glmer(formula=formula(paste(phenotype[p],FIX.formula,"+(",substr(RANDOM.formula,2,nchar(RANDOM.formula)),")")),data=z,family=eval(parse(text =family))))$coefficients)})
				betas<-lapply(time,function(w){lapply(model.coefs,function(z){(z)[which(rownames(z) %in% w ),]})})
				betas2<-matrix(unlist(lapply(betas,function(z)lapply(z,function(w){if(length(w)==0){w=NA}else{w}}))),ncol=length(time))
				betas3=cbind(betas2,betas3)
			}		
		}	
rownames(betas3)=names(model.coefs)
		

		#Step 1 for related subjects - Fitting Gene expressions:
		Pedigree.wide<-LongData$pedigree[!duplicated(LongData$ID)]
		model.coefs2=lapply(data.frame((EXPR)),function(z){as.data.frame(summary(lme(fixed=z~1,random=~1|Pedigree.wide))$coefficients$random)})	
		model.coefs<-do.call(cbind,model.coefs2)
		EXPR2<-(model.coefs)
		colnames(EXPR2)=names(model.coefs2)
		EXPR2<-t(EXPR2[order(rownames(EXPR2)),])	
		betas3<-betas3[order(rownames(betas3)),]

	}


	#Step 2: Analyses of Between Subject Variations
	GS=GS[order(rownames(GS)),]
	EXPR2=EXPR2[,order(colnames(EXPR2))]
	GS2= as.data.frame((GS))
	GS.sizes <- sapply(GS2,sum) # size of each gene set
	GS.data<-NULL	
	GS.data <- lapply(data.frame(GS2), function(z) (EXPR2)[z==1, ]); # creat data of each GS (rows=gene trends, columns=samples/families)
	GS.data <- lapply(GS.data,function(z) t(z)); # transfer genes in each GS (columns=gene trends, rows=samples/families)
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


      Cov.Pooled<-cov.shrink(betas3,verbose=FALSE, lambda.var=0);
      EIGEN.decom<-eigen(Cov.Pooled);  # eigen decomposition of pooled covariance for beta

	D<-EIGEN.decom$values;         
      U<-EIGEN.decom$vectors;
      cl<-t(betas3%*%U)/sqrt(D)          # adjust beta (rows=respons, columns=samples)


	# (3) T-like stats obtained on 'true' data
      sam.sumsquareT.obs  <- sapply(GS.data, function(z) T2.like.SAMGS(z,cl))
                                      # the T-like statistics obtained on 'true' data


     # (4) stats obtained on 'permuted' data
     sam.sumsquareT.permut <- matrix(NA,nbPermutations,nb.GeneSets)
     for(i in 1:nbPermutations) {

		if (familybased==FALSE){
      		ind <- sample(nb.Samples)
		}else{
			ind <- sample(nb.Families)
		}	
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
