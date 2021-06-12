# Propensity Score-Matching Methods for Nonexperimental Causal Studies 

setwd("~/workdir/statistics/Today-I-Learned/2021-1-AdvancedLinearModels")

library(MatchIt)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(gridExtra)
library(party)

data("lalonde")
table(lalonde$treat)
lalonde$race <- ifelse(lalonde$race == "black", 1, 0)
View(lalonde)

# 1. Matching with Logistic Regression 
# ==================================================================
# 1-1. difference in means : outcome variables 
lalonde %>%
  group_by(treat) %>%
  summarise(n_people = n(),
            mean_outcome = mean(re78),
            std_error = sd(re78) / sqrt(n_people))

with(lalonde, t.test(re78 ~ treat))


# 1-2. difference in means : pre-treatment covariates 
lalonde_cov <- c('age', 'educ', 'race', 'married', 'nodegree', 're74', 're75')
lalonde %>%
  group_by(treat) %>%
  select(one_of(lalonde_cov)) %>%
  summarise_all(funs(mean(., na.rm = T)))


# 2. Propensity Score Estimation 
fit1 <- glm(treat ~ age + educ + race + married + nodegree + 
              re74 + re75, data=lalonde, family=binomial())
lalonde$ps <- predict(fit1, lalonde, type="response")
summary(lalonde$ps[lalonde$treat==0])
summary(lalonde$ps[lalonde$treat==1])
summary(fit1)

# 2-2. Propensity Score Graph 
lalonde %>%
  mutate(ps.grp = round(ps/0.05) * 0.05) %>%
  group_by(treat, ps.grp) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  mutate(n2 = ifelse(treat == 0, yes = n, no = -1*n)) %>%
  ggplot(aes(x = ps.grp, y = n2, fill = as.factor(treat))) +
  geom_bar(stat = 'identity', position = 'identity') +
  geom_text(aes(label = n, x = ps.grp, y = n2 + ifelse(treat == 0, 8, -8))) +
  xlab('Probability of Job Training Program') +
  ylab(' ') +
  ggtitle('Propensity Score Distribution by Treatment Group') +
  scale_fill_discrete('') +
  scale_x_continuous(breaks = seq(0, 1, 0.05)) +
  theme(legend.position = 'bottom', legend.direction = 'vertical',
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())


# # Propensity Score Matching Graph 
# df1 <- lalonde[which(lalonde$treat==0),] ; df2 <- lalonde[which(lalonde$treat==1),] ; 
# df1 <- df1[order(df1$ps),] ; df2 <- df2[order(df2$ps),]
# rownames(df1) <- NULL ; rownames(df2) <- NULL
# 
# ggplot() + 
#   geom_line(data=df1, aes(x = as.numeric(row.names(df1)), y=ps), color='green') + 
#   geom_line(data=df2, aes(x = as.numeric(row.names(df2)), y=ps), color='red')
# 

# logistic regression, with all covariates 
m.out1 <- matchit(treat ~ age + educ + race + married 
                  + nodegree + re74 + re75, data = lalonde,
                  distance = "glm", link = "linear.probit")

m.out1
summary(m.out1)
round(summary(m.out1)$sum.all[,1:2],2)
plot(m.out1)
plot(summary(m.out1))
dim(match.data(m.out1)) # 370, 13 

plot(m.out1, type="jitter", interactive=FALSE)

dta.m1 <- match.data(m.out1)
View(dta.m1)


# 4.1 Visual inspection
fn_bal <- function(dta, variable) {
  dta$variable <- dta[, variable]
  dta$treat <- as.factor(dta$treat)
  support <- c(min(dta$variable), max(dta$variable))
  ggplot(dta, aes(x = distance, y = variable, color = treat)) +
    geom_point(alpha = 0.2, size = 1.3) +
    geom_smooth(method = "loess", se = F) +
    xlab("Propensity score") +
    ylab(variable) +
    theme_bw() +
    ylim(support)
}

grid.arrange(
  fn_bal(dta.m1, "age"),
  fn_bal(dta.m1, "educ") + theme(legend.position = "none"),
  fn_bal(dta.m1, "race"),
  fn_bal(dta.m1, "married") + theme(legend.position = "none"),
  fn_bal(dta.m1, "nodegree"),
  nrow = 3, widths = c(1, 0.8)
)


# 4.2 Difference in means 
dta.m1 %>%
  group_by(treat) %>%
  summarise(n_people = n(),
            mean_outcome = mean(re78),
            std_error = sd(re78) / sqrt(n_people))

with(dta.m1, t.test(re78 ~ treat))
lm_treat1 <- lm(re78 ~ treat, data = dta.m1)
summary(lm_treat1)

# 4.3 cov means 
dta.m1 %>%
  group_by(treat) %>%
  select(one_of(lalonde_cov)) %>%
  summarise_all(funs(mean))

lm_treat2 <- lm(re78 ~ treat + age + educ + race + married + nodegree + 
                  re74 + re75, data = dta.m1)
summary(lm_treat2)


# 2. Other Examples  
# =================================================================
data("lalonde")
table(lalonde$treat)
lalonde$race <- ifelse(lalonde$race == "black", 1, 0)


# ==================================================================
# (1) Caliper
m.out2 <- matchit(treat ~ age + educ + race + married + nodegree 
                  + re74 + re75, data = lalonde,
                  distance = "glm", caliper = .01)

m.out2
summary(m.out2)
round(summary(m.out2)$sum.all[,1:2],2)
plot(m.out2)
dim(match.data(m.out2)) # 180, 12 

dta.m2 <- match.data(m.out2)
plot(m.out2, type="jitter", interactive=FALSE)

grid.arrange(
  fn_bal(dta.m2, "age"), 
  fn_bal(dta.m2, "educ") + theme(legend.position = "none"),
  fn_bal(dta.m2, "race"),
  fn_bal(dta.m2, "married") + theme(legend.position = "none"),
  fn_bal(dta.m2, "nodegree"),
  nrow = 3, widths = c(1, 0.8)
)

with(dta.m2, t.test(re78 ~ treat))

lm_treat1 <- lm(re78 ~ treat, data = dta.m2)
summary(lm_treat1)

dta.m2 %>%
  group_by(treat) %>%
  select(one_of(lalonde_cov)) %>%
  summarise_all(funs(mean))

lm_treat2 <- lm(re78 ~ treat + age + educ + race + married + nodegree + 
                  re74 + re75, data = dta.m2)
summary(lm_treat2)


# ==================================================================
# (2) With Replacement
data("lalonde")
table(lalonde$treat)
lalonde$race <- ifelse(lalonde$race == "black", 1, 0)

m.out3 <- matchit(treat ~ age + educ + race + married 
                  + nodegree + re74 + re75, data = lalonde,
                  distance = "glm", replace = TRUE)

m.out3
summary(m.out3)
round(summary(m.out3)$sum.all[,1:2],2)
plot(m.out3)
dim(match.data(m.out3)) # 273, 13 
plot(m.out3, type="jitter", interactive=FALSE)


dta.m3 <- match.data(m.out3)

grid.arrange(
  fn_bal(dta.m3, "age"), 
  fn_bal(dta.m3, "educ") + theme(legend.position = "none"),
  fn_bal(dta.m3, "race"),
  fn_bal(dta.m3, "married") + theme(legend.position = "none"),
  fn_bal(dta.m3, "nodegree"),
  nrow = 3, widths = c(1, 0.8)
)

with(dta.m3, t.test(re78 ~ treat))

lm_treat1 <- lm(re78 ~ treat, data = dta.m3)
summary(lm_treat1)

dta.m3 %>%
  group_by(treat) %>%
  select(one_of(lalonde_cov)) %>%
  summarise_all(funs(mean))

lm_treat2 <- lm(re78 ~ treat + age + educ + race + married + nodegree + 
                  re74 + re75, data = dta.m3)
summary(lm_treat2)


# ==================================================================
# (3) 5:1 matching w/ replacement 
data("lalonde")
table(lalonde$treat)
lalonde$race <- ifelse(lalonde$race == "black", 1, 0)

m.out4 <- matchit(treat ~ age + educ + race + married 
                  + nodegree + re74 + re75, data = lalonde,
                  distance = "glm", replace = TRUE, ratio=5)

m.out4
summary(m.out4)
round(summary(m.out4)$sum.all[,1:2],2)
plot(m.out4)
dim(match.data(m.out4)) # 373, 13 

dta.m4 <- match.data(m.out4)
plot(m.out4, type="jitter", interactive=FALSE)

grid.arrange(
  fn_bal(dta.m4, "age"), 
  fn_bal(dta.m4, "educ") + theme(legend.position = "none"),
  fn_bal(dta.m4, "race"),
  fn_bal(dta.m4, "married") + theme(legend.position = "none"),
  fn_bal(dta.m4, "nodegree"),
  nrow = 3, widths = c(1, 0.8)
)

with(dta.m4, t.test(re78 ~ treat))

lm_treat1 <- lm(re78 ~ treat, data = dta.m4)
summary(lm_treat1)

dta.m4 %>%
  group_by(treat) %>%
  select(one_of(lalonde_cov)) %>%
  summarise_all(funs(mean))

lm_treat2 <- lm(re78 ~ treat + age + educ + race + married + nodegree + 
                  re74 + re75, data = dta.m4)
summary(lm_treat2)


# ==================================================================
# (4) Random Forest 
data("lalonde")
table(lalonde$treat)
lalonde$race <- ifelse(lalonde$race == "black", 1, 0)

library(randomForest)
m.out5 <- matchit(treat ~ age + educ + race + married 
                  + nodegree + re74 + re75, data = lalonde,
                  distance = "randomforest", replace = TRUE)

m.out5
summary(m.out5)
round(summary(m.out5)$sum.all[,1:3],2)
plot(summary(m.out5))
dim(match.data(m.out5)) # 253, 13 

dta.m5 <- match.data(m.out5)
plot(m.out5, type="jitter", interactive=FALSE)

grid.arrange(
  fn_bal(dta.m5, "age"), 
  fn_bal(dta.m5, "educ") + theme(legend.position = "none"),
  fn_bal(dta.m5, "race"),
  fn_bal(dta.m5, "married") + theme(legend.position = "none"),
  fn_bal(dta.m5, "nodegree"),
  nrow = 3, widths = c(1, 0.8)
)

with(dta.m5, t.test(re78 ~ treat))

lm_treat1 <- lm(re78 ~ treat, data = dta.m5)
summary(lm_treat1)

dta.m5 %>%
  group_by(treat) %>%
  select(one_of(lalonde_cov)) %>%
  summarise_all(funs(mean))

lm_treat2 <- lm(re78 ~ treat + age + educ + race + married + nodegree + 
                  re74 + re75, data = dta.m5)
summary(lm_treat2)

