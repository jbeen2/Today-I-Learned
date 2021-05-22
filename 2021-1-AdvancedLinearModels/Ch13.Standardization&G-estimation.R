library(readr)
library(dplyr)
library(tableone)
library(geepack)

# ======= Standardization ======= 
# \Sigma E(Y|L=l, A=a) * P(L=l)
# data : NHEFS from uncensored individuals (C=0)

fit <- glm(wt82_71 ~ qsmk + sex + race + age + I(age*age) + as.factor(education)
           + smokeintensity + I(smokeintensity*smokeintensity) + smokeyrs
           + I(smokeyrs*smokeyrs) + as.factor(exercise) + as.factor(active)
           + wt71 + I(wt71*wt71) + qsmk*smokeintensity, data=nhefs.nmw)
summary(fit)$coef

nhefs.nmw$predicted.meanY <- predict(fit, nhefs.nmw)  # y hat 
summary(nhefs.nmw$predicted.meanY[nhefs.nmw$cens==0]) # Min -10.88, Max 9.88 / Mean 2.64 


# create a dataset with 3 copies of each subject
nhefs.nmw$interv <- -1 # 1st copy: equal to original one

interv0 <- nhefs.nmw   # 2nd copy: treatment set to 0, outcome to missing
interv0$interv <- 0    # treatment A = 0 으로 고정 
interv0$qsmk <- 0
interv0$wt82_71 <- NA

interv1 <- nhefs.nmw   # 3rd copy: treatment set to 1, outcome to missing
interv1$interv <- 1    # treatment A = 1 으로 고정 
interv1$qsmk <- 1
interv1$wt82_71 <- NA

onesample <- rbind(nhefs.nmw, interv0, interv1) # combining datasets : 1629*3 = 4887 


# linear model to estimate mean outcome conditional on treatment and confounders
# parameters are estimated using original observations only (nhefs)
# parameter estimates are used to predict mean outcome for observations with
# treatment set to 0 (interv=0) and to 1 (interv=1)

std <- glm(wt82_71 ~ qsmk + sex + race + age + I(age*age)
           + as.factor(education) + smokeintensity
           + I(smokeintensity*smokeintensity) + smokeyrs
           + I(smokeyrs*smokeyrs) + as.factor(exercise)
           + as.factor(active) + wt71 + I(wt71*wt71) + I(qsmk*smokeintensity),
           data=onesample) # 3개의 dataset 으로 modeling 
summary(std) # coef 

onesample$predicted_meanY <- predict(std, onesample)

# estimate mean outcome in each of the groups interv=0, and interv=1
# this mean outcome is a weighted average of the mean outcomes in each combination
# of values of treatment and confounders, that is, the standardized outcome

mean(onesample[which(onesample$interv==-1),]$predicted_meanY) # 2.56

EY.a0.c0<-mean(onesample[which(onesample$interv==0),]$predicted_meanY)
EY.a1.c0<-mean(onesample[which(onesample$interv==1),]$predicted_meanY)
EY.a1.c0 - EY.a0.c0 # 3.52 



# ======= G-estimation =======
# estimation of denominator of ip weights for C
cw.denom <- glm(cens==0 ~ qsmk + sex + race + age + I(age^2)
                + as.factor(education) + smokeintensity + I(smokeintensity^2)
                + smokeyrs + I(smokeyrs^2) + as.factor(exercise)
                + as.factor(active) + wt71 + I(wt71^2),
                data = nhefs.nmw, family = binomial("logit"))
summary(cw.denom)

nhefs.nmw$pd.c <- predict(cw.denom, nhefs.nmw, type="response")
summary(nhefs.nmw$pd.c)
## Min. 1st Qu. Median Mean 3rd Qu. Max.
## 0.128 0.958 0.977 0.961 0.987 0.998

nhefs.nmw$wc <- ifelse(nhefs.nmw$cens==0, 1/nhefs.nmw$pd.c, NA) # observations with cens=1 only contribute to censoring models

# logit P(A=1|H(psi),L) = a0 + a1H(psi) + a2L
# Checking one possible value of psi
nhefs.nmw$psi <- 3.446  # 설명변수 처럼 넣어서, logistic 
nhefs.nmw$Hpsi <- nhefs.nmw$wt82_71 - nhefs.nmw$psi*nhefs.nmw$qsmk  # H(psi) = Y - psi * A
fit <- geeglm(qsmk ~ sex + race + age + I(age*age) + as.factor(education)
              + smokeintensity + I(smokeintensity*smokeintensity) + smokeyrs
              + I(smokeyrs*smokeyrs) + as.factor(exercise) + as.factor(active)
              + wt71 + I(wt71*wt71) + Hpsi, family=binomial, data=nhefs.nmw,
              weights=wc, id=seqn, corstr="independence")
tail(summary(fit)$coef)

# grid search for psi 
grid <- seq(from = 2,to = 5, by = 0.1)
j = 0
Hpsi.coefs <- cbind(rep(NA,length(grid)), rep(NA, length(grid)))
colnames(Hpsi.coefs) <- c("Estimate", "p-value")
for (i in grid){
  psi = i
  j = j+1
  nhefs.nmw$Hpsi <- nhefs.nmw$wt82_71 - psi * nhefs.nmw$qsmk
  gest.fit <- geeglm(qsmk ~ sex + race + age + I(age*age) + as.factor(education)
                     + smokeintensity + I(smokeintensity*smokeintensity) + smokeyrs
                     + I(smokeyrs*smokeyrs) + as.factor(exercise) + as.factor(active)
                     + wt71 + I(wt71*wt71) + Hpsi, family=binomial, data=nhefs.nmw,
                     weights=wc, id=seqn, corstr="independence")
  Hpsi.coefs[j,1] <- psi
  Hpsi.coefs[j,2] <- summary(gest.fit)$coefficients["Hpsi", "Pr(>|W|)"]
}

Hpsi.coefs

plot(Hpsi.coefs[,1], Hpsi.coefs[,2], type="b", xlab="psi", ylab="p-value",main="confidence interval of psi")
abline(h=0.05,lty="dashed")


##########################################################
# G-estimation: Closed form estimator linear mean models #
##########################################################
logit.est <- glm(qsmk ~ sex + race + age + I(age^2) + as.factor(education)
                 + smokeintensity + I(smokeintensity^2) + smokeyrs
                 + I(smokeyrs^2) + as.factor(exercise) + as.factor(active)
                 + wt71 + I(wt71^2), data = nhefs.nmw, weight = wc,
                 family = binomial())
summary(logit.est)

nhefs.nmw$pqsmk <- predict(logit.est, nhefs.nmw, type = "response") # coef 기반으로 확률값 
summary(nhefs.nmw$pqsmk)
## Min. 1st Qu. Median Mean 3rd Qu. Max.
## 0.051 0.178 0.243 0.262 0.325 0.789

# solve sum(w_c * H(psi) * (qsmk - E[qsmk | L])) = 0
# for a single psi and H(psi) = wt82_71 - psi * qsmk
# this can be solved as psi = sum( w_c * wt82_71 * (qsmk - pqsmk)) / sum(w_c * qsmk * (qsmk - pqsmk))
nhefs.c <- nhefs.nmw[which(!is.na(nhefs.nmw$wt82)),]
with(nhefs.c, sum(wc*wt82_71*(qsmk-pqsmk)) / sum(wc*qsmk*(qsmk - pqsmk)))  # 3.45 

#############################################################
# G-estimation: Closed form estimator for 2-parameter model #
#############################################################
diff = with(nhefs.c, qsmk - pqsmk)
diff2 = with(nhefs.c, wc * diff)
lhs = matrix(0,2,2)
lhs[1,1] = with(nhefs.c, sum(qsmk * diff2))
lhs[1,2] = with(nhefs.c, sum(qsmk * smokeintensity * diff2))
lhs[2,1] = with(nhefs.c, sum(qsmk * smokeintensity * diff2))
lhs[2,2] = with(nhefs.c, sum(qsmk * smokeintensity * smokeintensity * diff2))
rhs = matrix(0,2,1)
rhs[1] = with(nhefs.c, sum(wt82_71 * diff2))
rhs[2] = with(nhefs.c, sum(wt82_71 * smokeintensity * diff2))
psi = t(solve(lhs,rhs))
psi  # 2.86, 0.03 
