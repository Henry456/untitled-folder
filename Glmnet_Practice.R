library(glmnet)
data(QuickStartExample)
fit <- glmnet(x, y)
#dim(x)-->100 20. dim(y)-->100 1
plot(fit)
#axis above means that the number of nonzero coef
#at the current lambda

print(fit)#the function is to summary of glmnet path
coef(fit,s=0.1)

nx = matrix(rnorm(10*20),10,20)
predict(fit,newx=nx,s=c(0.1,0.05))

cvfit = cv.glmnet(x, y)
plot(cvfit)
coef(cvfit, s = "lambda.min") # or "lambda.1se"
predict(cvfit, newx = x[1:5,], s = "lambda.min")


fit = glmnet(x, y, alpha = 0.2, weights = c(rep(1,50),rep(2,50)), nlambda = 20)
print(fit)
plot(fit, xvar = "lambda", label = TRUE)


plot(fit, xvar = "dev", label = TRUE)
