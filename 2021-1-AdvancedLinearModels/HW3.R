library(readr)
library(dplyr)
library(tableone)
library(geepack)

nhefs <- read_csv("nhefs.csv")
View(nhefs)

# 82년도에 관측치가 없어 절단된 데이터 : censoring 
table(nhefs$qsmk, is.na(nhefs["wt82"]))
nhefs$cens <- ifelse(is.na(nhefs["wt82"]), 1, 0)

# Doubly Robust Estimator 

# 1. IPW for treatment A 
# Estimation of P(A=1|L) (=IP weights) : via logistic regression 

# A : qsmk (treatment)
fit <- glm(qsmk ~ sex + race + age + I(age^2) +
           as.factor(education) + smokeintensity +
           I(smokeintensity^2) + smokeyrs + I(smokeyrs^2) +
           as.factor(exercise) + as.factor(active) + wt71 + I(wt71^2),
           family = binomial(), data = nhefs)
fit  # coef 


# 2. Fit Outcome Regression 
# qsmk = 0 이면 1-p / qsmk=1 이면 p 
p.qsmk.obs <- ifelse(nhefs$qsmk == 0, 1-predict(fit, type = "response"),
                     predict(fit, type = "response"))

nhefs$w <- 1/p.qsmk.obs  # 1629 개의 weight : 1/{P(Ai=1|Li)}, 1/{P(Ai=0|Li)}
summary(nhefs$w)  # 1.056 ~ 16 

nhefs$r <- ifelse(nhefs$qsmk == 1, 1/p.qsmk.obs, -(1/p.qsmk.obs))
summary(nhefs$r)  # -2.92 ~ 16.01

# GEE, Y ~ A+R, weight <- IPW weight 
msm.w <- geeglm(wt82_71 ~ qsmk + r, data=nhefs, id=seqn, weights = w,  
                corstr="independence")
summary(msm.w)
# 
# # Use the predicted values from the outcome model
# nhefs$wt82_71_hat <- predict(msm.w, nhefs)
# 
# 
# # 3. Standardization
# # \Sigma E(Y|L=l, A=a) * P(L=l)
# # data : NHEFS from uncensored individuals (C=0)
# 
# y.fit <- glm(wt82_71_hat ~ qsmk + sex + race + age + I(age*age) + as.factor(education)
#            + smokeintensity + I(smokeintensity*smokeintensity) + smokeyrs
#            + I(smokeyrs*smokeyrs) + as.factor(exercise) + as.factor(active)
#            + wt71 + I(wt71*wt71) + qsmk*smokeintensity, data=nhefs)
# summary(y.fit)$coef

nhefs$predicted.meanY <- predict(msm.w, nhefs)  # y hat 
summary(nhefs$predicted.meanY[nhefs$cens==0]) # Min 1.07, Max 10.35 / Mean 2.58


# create a dataset with 3 copies of each subject
nhefs$interv <- -1 # 1st copy: equal to original one

interv0 <- nhefs   # 2nd copy: treatment set to 0, outcome to missing
interv0$interv <- 0    # treatment A = 0 으로 고정 
interv0$qsmk <- 0
interv0$predicted.meanY <- NA

interv1 <- nhefs   # 3rd copy: treatment set to 1, outcome to missing
interv1$interv <- 1    # treatment A = 1 으로 고정 
interv1$qsmk <- 1
interv1$predicted.meanY <- NA

onesample <- rbind(nhefs, interv0, interv1) # combining datasets : 1629*3 = 4887 


std <- glm(predicted.meanY ~ qsmk + sex + race + age + I(age*age)
           + as.factor(education) + smokeintensity
           + I(smokeintensity*smokeintensity) + smokeyrs
           + I(smokeyrs*smokeyrs) + as.factor(exercise)
           + as.factor(active) + wt71 + I(wt71*wt71) + I(qsmk*smokeintensity),
           data=onesample) # 3개의 dataset 으로 modeling 
summary(std) # coef 

onesample$predicted_meanY <- predict(std, onesample)
mean(onesample[which(onesample$interv==-1),]$predicted_meanY) 

EY.a0.c0<-mean(onesample[which(onesample$interv==0),]$predicted_meanY)
EY.a1.c0<-mean(onesample[which(onesample$interv==1),]$predicted_meanY)
EY.a1.c0 - EY.a0.c0   # 3.52 
