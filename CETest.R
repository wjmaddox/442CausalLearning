library(pcalg)

skcm <- read.csv("SKCM.csv",stringsAsFactors=FALSE)
skcm <- skcm[-which(skcm[,1]=="NO"),]
skcm <- skcm[!is.na(skcm[,1]),]

skcm2 <- skcm
for(x in 1:dim(skcm2)[2]){
  nas <- grep("\\[Not Applicable\\]",skcm[,x])
  nas2 <- grep("\\[Unknown\\]",skcm[,x])
  #nas3 <- grep("\\[Completed\\]",skcm[,x])
  skcm2[c(nas,nas2),x]<-NA
}

skcm.try <- data.frame(skcm2[,c(5,6,7,8,9,10)])

skcm.try[,1] <- as.numeric(as.factor(skcm.try[,1]))   
skcm.try[,3] <- as.numeric(as.factor(skcm.try[,3]))
skcm.try[,6] <- as.numeric(as.factor(skcm.try[,6]))
nas3 <- grep("\\[Completed\\]",skcm.try[,2])
skcm.try [nas3,2] <- NA
c <- cor(skcm.try,use="na.or.complete")

skcm.try[,2] <- as.numeric(skcm.try[,2])
skcm.try[,4] <- as.numeric(skcm.try[,4])
skcm.try[,5] <- as.numeric(skcm.try[,5])

skcm.fit <- pc(suffStat=list(C=c,n=nrow(skcm.try)),indepTest=gaussCItest,labels=colnames(skcm.try),alpha=.01)
iplotPC(skcm.fit)

skcm.fit2 <- rfci(suffStat=list(C=c,n=nrow(skcm.try)),indepTest=gaussCItest,labels=colnames(skcm.try),
                  alpha=.01,skel.method="stable.fast",conservative=TRUE)
iplotPC(skcm.fit2)

ida(3,2,cov(skcm.try,use="na.or.complete"),skcm.fit2,method="global",verbose=TRUE)
