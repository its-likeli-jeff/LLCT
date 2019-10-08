T2.like.SAMGS <-
function(DATA, cl){
  # DATA : expression data
  #     -> dataframe with  rows=genes,
  #                        columns=samples,
  # weight: weights to T2 statistics of genes. 
  # cl : response vectors for the samples IN THE SAME ORDER AS IN DATA
  #    -> a dataframe with rows=respons,
  #                      colmns=samples
 
   cl<-as.matrix(cl)
   DATA<-as.matrix(DATA)
   Sigma<-cov(t(DATA),t(cl))
   dim.pq<-dim(Sigma)
   if(dim.pq[1]>dim.pq[2]) max(eigen(t(Sigma)%*%Sigma,symmetric=T,only.values=T)[[1]])
   else max(eigen(Sigma%*%t(Sigma),symmetric=T,only.values=T)[[1]])
}
