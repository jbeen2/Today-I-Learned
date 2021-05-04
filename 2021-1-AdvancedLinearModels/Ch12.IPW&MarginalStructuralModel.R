# ===============  Goal =============== 
# To estimate the ACE of smoking cessation (treatment) A 
# on on weight gain (outcome) Y 

library(readr)
library(dplyr)
library(tableone)
library(geepack)

nhefs.nmw <- read_csv("nhefs.csv")
lm(wt82_71 ~ qsmk, data = nhefs.nmw)  # E(Y|A=0) = 1.984 / E(Y|A=1) = 1.984+2.541

# smoking cessation 
predict(lm(wt82_71 ~ qsmk, data = nhefs.nmw), data.frame(qsmk=1))  # 4.5250

# No smoking cessation 
predict(lm(wt82_71 ~ qsmk, data = nhefs.nmw), data.frame(qsmk=0))  # 1.9845 

vars<-c("age","wt71","smokeintensity","smokeyrs","sex","race","education","exercise","active")
cat.var<-c("sex","race","education","exercise","active")  # categorical variables 
nhefs.nmw %>% mutate_at(cat.var, funs(factor(.)))         # convert cat.var to factor

# <CreateTableOne> 정리 잘 해 주는 함수 
# continuous variables : 평균/표준편차, t-test 
# discrete variables : confusion matrix, Chi-square test 
tab <- CreateTableOne(vars = vars, strata = "qsmk" ,
                      data = nhefs.nmw, factorVars = cat.var)
print(tab, showAllLevels = TRUE, formatOptions = list(big.mark = ","))   # Show All Levels! 


# Age : "Confounder" of the effect of A on Y 

# ========== PROCESS ==========
# 1. P(A=1|L) -> weight (via logistic regression)
# 2. GEE weight , Y ~ A 
# =============================

# IPW for treatment A 
# Estimation of P(A=1|L) (=IP weights) : via logistic regression 
fit <- glm(qsmk ~ sex + race + age + I(age^2) +
             as.factor(education) + smokeintensity +
             I(smokeintensity^2) + smokeyrs + I(smokeyrs^2) +
             as.factor(exercise) + as.factor(active) + wt71 + I(wt71^2),
             family = binomial(), data = nhefs.nmw)
fit

p.qsmk.obs <- ifelse(nhefs.nmw$qsmk == 0, 1 - predict(fit, type = "response"),
                     predict(fit, type = "response"))

# IPW
nhefs.nmw$w <- 1/p.qsmk.obs  # weight : 1/{P(Ai=1|Li)}, 1/{P(Ai=0|Li)}
summary(nhefs.nmw$w)  # 1.056 ~ 16 
# 16 값을 가진 사람은 다른 사람들보다 조금 남다르다 
# 한 사람의 data 가 16 배로 뻥튀기 되어서 pseudo population 에 들어간다 
# 이를 통해 두 집단의 비율이 맞춰진다! 

# GEE (Generalized Estimating Equations) 
# y들끼리의 correlation structure 틀리더라도, beta consistency 잘 구할 수 있다 
# with Sandwich Covariance Matrix 
msm.w <- geeglm(wt82_71 ~ qsmk, data=nhefs.nmw, weights=w, id=seqn,
                corstr="independence")
summary(msm.w) # E(Y|A=0) = 1.75 / E(Y|A=1) = 1.75+3.52 : difference causal effect 의 estimator 값 

beta <- coef(msm.w) ; SE <- coef(summary(msm.w))[,2]
lcl <- beta-qnorm(0.975)*SE ; ucl <- beta+qnorm(0.975)*SE
cbind(beta, lcl, ucl)   # 나이대를 맞춰 보았더니, 3.5 년 정도로 격차가 더 벌어졌다 


# ======= Stabilized IPW =======
# A 의 종류마다 다른 p를 주는 방법, group 별로 weight 다르게 부여  

# estimation of denominator of ip weights
denom.fit <- glm(qsmk ~ as.factor(sex) + as.factor(race) + age + I(age^2) +
                 as.factor(education) + smokeintensity +
                 I(smokeintensity^2) + smokeyrs + I(smokeyrs^2) +
                 as.factor(exercise) + as.factor(active) + wt71 + I(wt71^2),
                 family = binomial(), data = nhefs.nmw)
#denom.fit
pd.qsmk <- predict(denom.fit, type = "response")
summary(pd.qsmk)  # Min 0.053 ~ Max 0.793 

# estimation of numerator of ip weights
numer.fit <- glm(qsmk~1, family = binomial(), data = nhefs.nmw)
numer.fit  # -1.03 

pn.qsmk <- predict(numer.fit, type = "response")
nhefs.nmw$sw <- ifelse(nhefs.nmw$qsmk == 0, ((1-pn.qsmk)/(1-pd.qsmk)),
                       (pn.qsmk/pd.qsmk))
summary(nhefs.nmw$sw)  # Min 0.331 ~ Max 4.205 

msm.sw <- geeglm(wt82_71 ~ qsmk, data=nhefs.nmw, weights=sw, id=seqn,
                 corstr="independence")
summary(msm.sw)        # E(Y|A=0) = 1.75 / E(Y|A=1) = 1.75+3.52 

beta2 <- coef(msm.sw) ; SE2 <- coef(summary(msm.sw))[,2]
lcl2 <- beta2-qnorm(0.975)*SE2 ; ucl2 <- beta2+qnorm(0.975)*SE2
cbind(beta2, lcl2, ucl2)   # 결과가 좀 더 안정적 ... 



# ======= Marginal Structural Models =======
# E(Y|A) = theta0 + theta1 * A + theta2 * A^2 

# Analysis restricted to subjects reporting <=25 cig/day at baseline
nhefs.nmv.s <- subset(nhefs.nmw, smokeintensity <=25)

# estimation of denominator of ip weights
den.fit.obj <- lm(smkintensity82_71 ~ as.factor(sex) +
                    as.factor(race) + age + I(age^2) +
                    as.factor(education) + smokeintensity + I(smokeintensity^2) +
                    smokeyrs + I(smokeyrs^2) + as.factor(exercise) + as.factor(active) + wt71 +
                    I(wt71^2), data = nhefs.nmv.s)
p.den <- predict(den.fit.obj, type = "response")  # estimate prob 
dens.den <- dnorm(nhefs.nmv.s$smkintensity82_71, p.den, summary(den.fit.obj)$sigma) # normalize 

# estimation of numerator of ip weights
num.fit.obj <- lm(smkintensity82_71 ~ 1, data = nhefs.nmv.s)
p.num <- predict(num.fit.obj, type = "response")
dens.num <- dnorm(nhefs.nmv.s$smkintensity82_71, p.num, summary(num.fit.obj)$sigma)

nhefs.nmv.s$sw.a = dens.num/dens.den
summary(nhefs.nmv.s$sw.a)  # Min 0.19 ~ Max 4.98

# GEE 
msm.sw.cont <- geeglm(wt82_71 ~ smkintensity82_71 + I(smkintensity82_71*smkintensity82_71),
                      data=nhefs.nmv.s, weights=sw.a, id=seqn, corstr="independence")
summary(msm.sw.cont)  # 2.01 + (-0.11) * A + (0.0023) * A^2 

beta <- coef(msm.sw.cont) ; SE <- coef(summary(msm.sw.cont))[,2]
lcl <- beta-qnorm(0.975)*SE ; ucl <- beta+qnorm(0.975)*SE
cbind(beta, lcl, ucl)


# ======= Effect Modifier ======= 
# (Sex) A 에 영향을 주지 않지만 Y 에 영향을 주는 것 
table(nhefs.nmw$sex)

# estimation of denominator of ip weights
denom.fit <- glm(qsmk ~ as.factor(sex) + as.factor(race) + age + I(age^2) +
                   as.factor(education) + smokeintensity +
                   I(smokeintensity^2) + smokeyrs + I(smokeyrs^2) +
                   as.factor(exercise) + as.factor(active) + wt71 + I(wt71^2),
                 family = binomial(), data = nhefs.nmw)
summary(denom.fit)
pd.qsmk <- predict(denom.fit, type = "response")

# estimation of numerator of ip weights
numer.fit <- glm(qsmk~as.factor(sex), family = binomial(), data = nhefs.nmw)
summary(numer.fit)

pn.qsmk <- predict(numer.fit, type = "response")
nhefs.nmw$sw.a <- ifelse(nhefs.nmw$qsmk == 0, ((1-pn.qsmk)/(1-pd.qsmk)),
                         (pn.qsmk/pd.qsmk))  # SW = f(A|V) / f(A|V,L)
summary(nhefs.nmw$sw.a)  # Min 0.29 ~ Max 3.68 


# ======= Selection Bias =======
# Censoring : 중간에 도망간 그룹! 
# 모두가 관찰되었다고 가정했을 때, A treatment 의 차이를 보고 싶다 
table(nhefs.nmw$qsmk, nhefs.nmw$death)
summary(nhefs.nmw[which(nhefs.nmw$death==0),]$wt71)
summary(nhefs.nmw[which(nhefs.nmw$death==1),]$wt71)

# estimation of denominator of ip weights for C
denom.cens <- glm(death ~ as.factor(qsmk) + as.factor(sex) +
                    as.factor(race) + age + I(age^2) +
                    as.factor(education) + smokeintensity +
                    I(smokeintensity^2) + smokeyrs + I(smokeyrs^2) +
                    as.factor(exercise) + as.factor(active) + wt71 + I(wt71^2),
                    family = binomial(), data = nhefs.nmw)
pd.cens <- 1-predict(denom.cens, type = "response")

# estimation of numerator of ip weights for C
numer.cens <- glm(death~as.factor(qsmk), family = binomial(), data = nhefs.nmw)
summary(numer.cens)  # Wc = f(C|A) / f(C|A,L)
pn.cens <- 1-predict(numer.cens, type = "response")

nhefs.nmw$sw.a <- ifelse(nhefs.nmw$qsmk == 0, ((1-pn.qsmk)/(1-pd.qsmk)),
                     (pn.qsmk/pd.qsmk))
nhefs.nmw$sw.c <- pn.cens/pd.cens
nhefs.nmw$sw <- nhefs.nmw$sw.c*nhefs.nmw$sw.a
summary(nhefs.nmw$sw)  # Min 0.28 ~ Max 20.65 

msm.sw <- geeglm(wt82_71~qsmk, data=nhefs.nmw,
                 weights=sw, id=seqn, corstr="independence")
beta <- coef(msm.sw) ; SE <- coef(summary(msm.sw))[,2]
lcl <- beta-qnorm(0.975)*SE ; ucl <- beta+qnorm(0.975)*SE
cbind(beta, lcl, ucl)
