#################################################################33
#merge clinical, followup & nte event data
####################################################################
rAndF <- function(df){#makes first row the column names
  colnames(df) <- as.character(unname(unlist(df[1,])))
  df <- df[-(1:3),]
  return(df)
}

nte.loc <- "nationwidechildrens.org_clinical_nte_kirc.txt"
followup.loc <- "nationwidechildrens.org_clinical_follow_up_v1.0_kirc.txt"
clinical.loc <- "nationwidechildrens.org_clinical_patient_kirc.txt"

kirc.allNC <- read.delim(clinical.loc,header=FALSE,sep="\t",stringsAsFactors=FALSE)
kirc.followupNC <- read.delim(followup.loc,header=FALSE,sep="\t",stringsAsFactors=FALSE)
kirc.nteNC <- read.delim(nte.loc,header=FALSE,sep="\t",stringsAsFactors=FALSE)

kirc.all <- rAndF(kirc.allNC)
kirc.followup <- rAndF(kirc.followupNC)
kirc.NTE <- rAndF(kirc.nteNC)

followups <- which(as.character(kirc.all[1,])%in%as.character(kirc.followup[1,]))

kirc <-kirc.all
row.names(kirc) <- 1:dim(kirc)[1]
kirc$FollowUp <- "NO"

in.CF <- which(colnames(kirc.followup)%in%colnames(kirc))

for (p in 1:dim(kirc.followup)[1]){
  w <- which(as.character(kirc[,1])==as.character(kirc.followup[p,1]))
  if (length(w)!=0){
    #P <- c(P,p)
    cf.tmp <- kirc.followup[p,in.CF]
    cf.tmp$FollowUp <- "YES"
    for (c in 1:length(colnames(cf.tmp))){
      w2 <- which(colnames(kirc)==colnames(cf.tmp)[c])
      if(!is.na(cf.tmp[1,c])){
        kirc[w,w2]<-as.character(cf.tmp[1,c])
      }
    }
  }
  else{
    #nP <- c(nP,p)
    cf.tmp <-kirc.followup[p,in.CF]
    cf.tmp$FollowUp <- "NEW"
    kirc <- rbind.fill(kirc,cf.tmp)
  }
  #print(p)
}

kirc2 <- kirc
w3 <- which(as.character(kirc.NTE$bcr_patient_barcode)%in%as.character(kirc$bcr_patient_barcode))
nulldf <- data.frame(matrix(ncol=(ncol(kirc.NTE)-1),nrow=dim(kirc2)[1]))
colnames(nulldf)<-colnames(kirc.NTE)[-1]
kirc2 <- cbind(kirc2,nulldf)
for (p in 1:dim(kirc.NTE)[1]){
  NTE.id <- as.character(kirc.NTE[p,1])
  w4 <- which(as.character(kirc2$bcr_patient_barcode)==NTE.id)
  kirc2[w4,(ncol(kirc)+1):(ncol(kirc2))] <- kirc.NTE[p,-1]
}
