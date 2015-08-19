###########################################################3
#Name:R_ClnicalMergingMissing.R
#Date: 3/27
#Author: WM
#KIRC KIRP LGG LUAD LUSC PRAD SKCM STAD
############################################################
tumor2 <- c("KIRC", "KIRP", "LGG", "LUAD", "LUSC", "PRAD", "SKCM", "STAD")
for (tumor in tumor2){
  loc <- paste("/fs2/tom/ClinicalData/",toupper(tumor),"/Complete_Clinical_Set/NCH__Biotab/Level_2/nationwidechildrens.org_clinical",sep="")
  all.loc <- system(paste("ls ",loc,"_patient_",tolower(tumor),".txt",sep=""),intern=TRUE)#clinical_patient
  followupsTMP.loc <- system(paste("ls ",loc,"_follow_up*",sep=""),intern=TRUE)#both patient & nte followup
  
  nteFollowups.loc <- followupsTMP.loc[grep("nte",followupsTMP.loc)]
  nteFollowup.loc <- nteFollowups.loc[length(nteFollowups.loc)]
  
  followups.loc <- followupsTMP.loc[-grep("nte",followupsTMP.loc)]
  followup.loc <- followups.loc[length(followups.loc)]
  nte.loc <- system(paste("ls ",loc,"_nte*",sep=""),intern=TRUE)
  
  if(length(all.loc)!=0){
    lihc.allNC <- read.delim(all.loc,header=FALSE,sep="\t",stringsAsFactors=FALSE)
    lihc.all <- rAndF(lihc.allNC)
  }  else{lihc.all <- NULL}
  if(length(followup.loc)!=0){
    lihc.followupNC <- read.delim(followup.loc,header=FALSE,sep="\t",stringsAsFactors=FALSE)
    lihc.followup <- rAndF(lihc.followupNC)
  }  else{lihc.followup <- NULL}
  if(length(nte.loc)!=0){
    lihc.nteNC <- read.delim(nte.loc,header=FALSE,sep="\t",stringsAsFactors=FALSE)
    lihc.nte <- rAndF(lihc.nteNC)
  }  else{lihc.nte <- NULL}
  if(length(nteFollowup.loc)!=0){
    lihc.nteFollowupNC <- read.delim(nteFollowup.loc,header=FALSE,sep="\t",stringsAsFactors=FALSE)
    lihc.nteFollowup <- rAndF(lihc.nteFollowupNC)
  }  else{lihc.nteFollowup <- NULL}
  
  print("successfully read in files")
  ###############
  ##merge lihc.all & lihc.followup
  followups <- which(as.character(lihc.all[1,])%in%as.character(lihc.followup[1,]))
  
  lihc <-lihc.all
  row.names(lihc) <- 1:dim(lihc)[1]
  lihc$FollowUp <- "NO"
  
  in.CF <- which(colnames(lihc.followup)%in%colnames(lihc))
  #in.lihc <- which(colnames(lihc)%in%colnames(lihc.followup))
  
  #nP<-c()
  #P<-c()
  if(!is.null(lihc.followup)){
    for (p in 1:dim(lihc.followup)[1]){
      w <- which(as.character(lihc[,1])==as.character(lihc.followup[p,1]))
      if (length(w)!=0){
        #P <- c(P,p)
        cf.tmp <- lihc.followup[p,in.CF]
        cf.tmp$FollowUp <- "YES"
        for (c in 1:length(colnames(cf.tmp))){
          w2 <- which(colnames(lihc)==colnames(cf.tmp)[c])
          if(!is.na(cf.tmp[1,c])){
            lihc[w,w2]<-as.character(cf.tmp[1,c])
          }
        }
      }
      else{
        #nP <- c(nP,p)
        cf.tmp <-lihc.followup[p,in.CF]
        cf.tmp$FollowUp <- "NEW"
        lihc <- rbind.fill(lihc,cf.tmp)
      }
      #print(p)
    }
  }
  ####################
  #merged nte & nte followup
  
  lihc.NTE <-lihc.nte
  lihc.NTE$FollowUpNTE <- "NO"
  
  in.CFNTE <- which(colnames(lihc.nteFollowup)%in%colnames(lihc.NTE))
  #in.clinical <- which(colnames(clinical)%in%colnames(c.F))
  if(!(is.null(dim(lihc.nteFollowup)))){
    for (p in 1:dim(lihc.nteFollowup)[1]){
      w <- which(as.character(lihc.NTE[,1])==as.character(lihc.nteFollowup[p,1]))
      if (length(w)!=0){
        cf.tmp <- lihc.nteFollowup[p,in.CFNTE]
        cf.tmp$FollowUpNTE <- "YES"
        for (c in 1:length(colnames(cf.tmp))){
          w2 <- which(colnames(lihc.NTE)==colnames(cf.tmp)[c])
          if(!is.na(cf.tmp[1,c])){
            lihc.NTE[w,w2]<-cf.tmp[1,c]
          }
        }
      }
      else{
        cf.tmp <-lihc.nteFollowup[p,in.CFNTE]
        cf.tmp$FollowUpNTE <- "NEW"
        lihc.NTE <- rbind.fill(lihc.NTE,cf.tmp)
      }
      #print(p)
    }
  }
  lihc2 <- lihc
  if(!is.null(dim(lihc.NTE))){
    w3 <- which(as.character(lihc.NTE$bcr_patient_barcode)%in%as.character(lihc$bcr_patient_barcode))
    nulldf <- data.frame(matrix(ncol=(ncol(lihc.NTE)-1),nrow=dim(lihc2)[1]))
    colnames(nulldf)<-colnames(lihc.NTE)[-1]
    lihc2 <- cbind(lihc2,nulldf)
    for (p in 1:dim(lihc.NTE)[1]){
      NTE.id <- as.character(lihc.NTE[p,1])
      w4 <- which(as.character(lihc2$bcr_patient_barcode)==NTE.id)
      lihc2[w4,(ncol(lihc)+1):(ncol(lihc2))] <- lihc.NTE[p,-1]
    }
  }
  for(x in 1:dim(lihc2)[2]){
    nas <- grep("\\[Not Available\\]",as.character(lihc2[,x]))
    lihc2[nas,x]<-NA
  }
  print(paste("Outputting File: ",toupper(tumor),".csv",sep=""))
  write.csv(lihc2,paste(toupper(tumor),".csv",sep=""),row.names=FALSE,quote=FALSE)
}