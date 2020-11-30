library("mediation")
set.seed(2020)
data("framing", package = "mediation")

# Bayesian 
med.fit <- lm(emo ~  treat + age + educ + gender + income, data = framing)
out.fit <- lm(cong_mesg ~ emo + treat + age + educ + gender + income, data = framing)
med.out <- mediate(med.fit, out.fit, treat="treat", mediator="emo")
summary(med.out)
plot(med.out)


# Nonparametric Bootstrap 
med.out2 <- mediate(med.fit, out.fit, boot = TRUE, treat="treat", mediator="emo")
summary(med.out2)
plot(med.out2)
