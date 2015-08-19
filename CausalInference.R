########################################
#build models and plot graphs
########################################

#########################
#read in data
load("graph.RData")
load("graph2.RData")

load("skcmSub.RData")
skcm.int$new_tumor_event_dx_indicator <- as.numeric(as.character(skcm.int$new_tumor_event_dx_indicator))-1
skcm.int$primary_multiple_at_dx <- as.numeric(as.character(skcm.int$primary_multiple_at_dx))-1
w <- which(is.na(skcm.int),arr.ind=T)
skcm.int <- skcm.int[-w[,1],]
#make binary data
tkplot(skcm.fi2G)
tkplot(int.fitG)

#####################################################
#marginal structural models for the simple stuff
#simple logistic models for probability of metastasis
#P(new_tumor_event_dx_indicator) ~ primary_multiple_at_dx
w1 <- which(skcm.int$primary_multiple_at_dx==1)
w2 <- which(skcm.int$primary_multiple_at_dx==0)
#model2.p <- glm(new_tumor_event_dx_indicator ~ primary_multiple_at_dx, data=skcm.int[w1,],family=binomial(),weight=sw[w1])
#model1.p <- glm(new_tumor_event_dx_indicator ~ primary_multiple_at_dx, data=skcm.int[w2,],family=binomial(),weight=sw[w2])

#covariates are rest of dataset
denom.fit <- glm(primary_multiple_at_dx ~ race+vital_status+tumor_status+retrospective_collection,family=binomial(),data=skcm.int)
denom.pT <- predict(denom.fit,type="response")

num.fit <- glm(primary_multiple_at_dx ~ 1, family=binomial(), data=skcm.int)
num.p <- predict(num.fit,type="response")

sw <- ifelse(skcm.int$primary_multiple_at_dx==0,((1-as.numeric(num.p))/(1-as.numeric(denom.pT))),as.numeric(num.p)/as.numeric(denom.pT))

#model1.p <- glm(new_tumor_event_dx_indicator ~ primary_multiple_at_dx, data=skcm.int,family=binomial(),weight=sw)
model.1 <- glm(new_tumor_event_dx_indicator ~ primary_multiple_at_dx, data=skcm.int,family=binomial(),weight=sw)

#primary multiple at diagnosis seems to be a perfect predictor of having a new tumor event...

#P(nte_dx) ~ race
ce.2 <- glm(new_tumor_event_dx_indicator ~ race, data=skcm.int, family=binomial())
#P(nte_dx) ~ race +primary_mult
denom.fit <- lm(race ~ primary_multiple_at_dx+vital_status+tumor_status+retrospective_collection,family=gaussian(),data=skcm.int)
denom.pT <- predict(denom.fit,type="response")
denom.dens <- dnorm(skcm.int$race,denom.pT,summary(denom.fit)$sigma)

num.fit <- lm(race ~ 1, family=gaussian(), data=skcm.int)
num.p <- predict(num.fit,type="response")
num.dens <- dnorm(skcm.int$race,num.p,summary(num.fit)$sigma)

#sw <- ifelse(skcm.int$primary_multiple_at_dx==0,((1-as.numeric(num.p))/(1-as.numeric(denom.pT))),as.numeric(num.p)/as.numeric(denom.pT))
sw <- num.dens/denom.dens

model.2 <- glm(new_tumor_event_dx_indicator ~ race, data=skcm.int, family=binomial(), weights=sw)
ce.race <- as.numeric(coef(model.2))[2]

##################################################
#instrumental variable analysis for tumor_status & retrospective collection

#vital_status is the instrument
#vital_status is binary, thus estimatior is simple
#(e(y|z=0)-e(y|z=1))/(e(x|z=0)-e(y|z=1))
z0 <- which(skcm.int$vital_status==1)#1 is alive
z1 <- which(skcm.int$vital_status==2)#2 is dead
ey.z0 <- mean(skcm.int$new_tumor_event_dx_indicator[z0])
ey.z1 <- mean(skcm.int$new_tumor_event_dx_indicator[z1])

ex1.z0 <- mean(skcm.int$tumor_status[z0])
ex1.z1 <- mean(skcm.int$tumor_status[z1])

ex2.z0 <- mean(skcm.int$retrospective_collection[z0])
ex2.z1 <- mean(skcm.int$retrospective_collection[z1])

ts.lvae <- (ey.z0-ey.z1)/(ex1.z0-ex1.z1)
rc.lvae <- (ey.z0-ey.z1)/(ex2.z0-ex2.z1)
