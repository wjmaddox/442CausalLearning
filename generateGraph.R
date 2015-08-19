#data processing for SKCM TCGA dataset
#SKCM data is cleaned and merged version of publicly available TCGA SKCM data

library(pcalg)
#pcalg is R package for generating causal graphs using the PC algorithm

#read in data
skcm <- read.csv("SKCM.csv",stringsAsFactors=FALSE)
skcm <- skcm[-which(skcm[,1]=="NO"),]
skcm <- skcm[!is.na(skcm[,1]),]

#ensure that missing data is ignored
skcm2 <- skcm
for(x in 1:dim(skcm2)[2]){
  nas <- grep("\\[Not Applicable\\]",skcm[,x])
  nas2 <- grep("\\[Unknown\\]",skcm[,x])
  #nas3 <- grep("\\[Completed\\]",skcm[,x])
  if(length(c(nas,nas2))>0){skcm2[c(nas,nas2),x]<-NA}
}

#remove columns of data that have mostly missing data
skcm2 <- skcm2[-c(1,2,3,4,57,56,36,38,45:48,49,50:54,60,68,71,59,57,35,25,58)]

#remove factors
for (x in 1:ncol(skcm2)){
  skcm2[,x] <- as.numeric(as.factor(skcm2[,x]))
}


#42: 2 is yes on metastasis
c <- cor(skcm2,use="pairwise.complete.obs")
c[is.na(c)]<-1e-10


skcm.fixedEdges <- matrix(0,nrow=ncol(skcm2),ncol=ncol(skcm2))
rownames(skcm.fixedEdges) <- colnames(skcm2)
colnames(skcm.fixedEdges) <- colnames(skcm2)

#manually create edges for obvious cases (ie metastasis location -> metastasis)
skcm.fixedEdges[13,11] <- 1
skcm.fixedEdges[11,13] <- 1
skcm.fixedEdges[c(35:36,38:49),34] <- 1
skcm.fixedEdges[34,c(35:36,38:49)] <- 1

#need to manually remove poor data values in the data frame
#add in the fixed edges
#create causal graph
skcm.fit <- pc(suffStat=list(C=c,n=nrow(skcm2)),fixedEdges=skcm.fixedEdges,indepTest=gaussCItest,labels=colnames(skcm2),alpha=.25)
#plot graph
iplotPC(skcm.fit@graph)

skcm.am <- showAmat(skcm.fit)
rs <- rowSums(skcm.am)
rs <- rs[-which(as.numeric(rs)==0)]

#run again - this time setting all correlations to non-zero
#this produces better results bc more data is allowed to be included
skcm3 <- skcm2[names(rs)]
c2 <- cor(skcm3,use="pairwise.complete.obs")
c2[is.na(c2)]<-1e-10
skcm.fi2 <- pc(suffStat=list(C=c2,n=nrow(skcm3)),fixedEdges=skcm.fixedEdges[names(rs),names(rs)],indepTest=gaussCItest,labels=colnames(skcm3),alpha=.05)
iplotPC(skcm.fi2@graph)

skcm.fi2G <- igraph.from.graphNEL(skcm.fi2@graph)
save(skcm.fi2G,file="graph.RData")
####################
# 
# #skcm.fit <- pc(suffStat=list(dm=skcm2))
# lf <- c()
# for (i in 1:ncol(skcm2)){
#   lf <- c(lf,length(levels(as.factor(skcm2[,i]))))
# }
# dm.ss <- list(dm=skcm2,nlev=lf,adaptDF=FALSE)
# #skcm.fit2 <- rfci(suffStat=list(C=c,n=nrow(skcm2)),indepTest=gaussCItest,labels=colnames(skcm2),
# #                  alpha=.01,skel.method="stable.fast",conservative=TRUE)
# skcm.fi2 <- pc(suffStat=dm.ss,indepTest=disCItest,labels=colnames(skcm2),alpha=.01)
# iplotPC(skcm.fit2)
# 
# ida(3,2,cov(skcm.try,use="pairwise.complete.obs"),skcm.fit2,method="global",verbose=TRUE)

#####################
#glm(new_tumor_event_dx_indicator ~ ., data=skcm3,family="binomial")
v <- 1:42[-27] #27 is new_tumor_event_dx
ida.ce <- list(NULL)
for (x in v){
  ida.ce[[x]] <- ida(x,27,graphEst=skcm.fi2@graph,mcov=cov(skcm3,use="pairwise.complete.obs"))
}
#now a list of possible causal effects
ida.ceMax <- lapply(ida.ce,max)
ida.ceMax <- unlist(ida.ceMax)
names(ida.ceMax) <- names(skcm3)[v]
ida.rank <- ida.ceMax[order(as.numeric(ida.ceMax))]
ida.interest <- names(tail(ida.rank))
print(ida.interest)

##########################################
#only take subgraph of these six interesting types
#highest possible causal effects on the metastasis probability
skcm.int <- skcm3[ida.interest]
c3 <- cor(skcm.int,use="pairwise.complete.obs")
int.fit <- pc(suffStat=list(C=c3,n=nrow(skcm.int)),indepTest=gaussCItest,labels=colnames(skcm.int),alpha=.05)
plot(int.fit)
int.fitG <- igraph.from.graphNEL(int.fit@graph)
save(int.fitG,file="graph2.RData")
save(skcm.int,file="skcmSub.RData")
