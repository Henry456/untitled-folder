library(penalized)
library(survival)
data(nki70)
set.seed(1)
fit <- penalized(ER, ~DIAPH3+NUSAP1, data=nki70, lambda1 = 1)
coefficients(fit, "all")
loglik(fit)
penalty(fit)

predict(fit, ~DIAPH3+NUSAP1, data=nki70[1:3,])
predict(fit, nki70[1:3, c("DIAPH3", "NUSAP1")])

pred <- predict(fit, nki70[1:3, c("DIAPH3", "NUSAP1")])
survival(pred, time=5)
coefficients(fit)
coefficients(fit, standardize = TRUE)
weights(fit)
fit <- penalized(Surv(time, event), nki70[, 8:77], ~ER, lambda2=1)
fit <- penalized(Surv(time, event)~ER, nki70[,8:77], lambda2=1)
fit <- penalized(Surv(time,event), nki70[,8:77], lambda1=1, steps=50, trace=FALSE)


fit <- cvl(Surv(time,event), nki70[,8:77], lambda1=1, fold=10)
fit$cvl

fit2 <- cvl(Surv(time,event), nki70[,8:77], lambda2=10, approximate=TRUE)
plot(fit1$lambda, fit1$cvl, type="l")
plot(fit2$lambda, fit2$cvl, type="l", log="x")

