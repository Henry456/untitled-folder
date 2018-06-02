#1)Generate required data set
#2)Compare Compare flare, penalized, and glmnet for sparse logistic regression
  #on above data
library(glmnet)
## generate data
n = 50
d = 100
set.seed(1)
X = matrix(rnorm(n*d), n, d)
beta = c(1, 1, 1, 1,rep(0,d-4))
Y = rep(0, 50)
#Generate the 
for(i in 1:50){
  Y[i] = exp(X[i, ] %*% beta)/(1+exp(X[i, ] %*% beta))
  if(Y[i] > 0.5) Y[i] = 1
  else Y[i] = 0
}

#apply the cv.glmnet method 
grid <- 10^seq(10,-2,length=100)
#This uses the cross-validation
cv.fit = cv.glmnet(X, Y, family = "binomial", lambda = grid)
#We can get the estimate for beta from coef()
lasso.coeff1 = coef(cv.fit)
#length(lasso.coeff1[lasso.coeff1 != 0])
#[1] 8
plot(cv.fit)
print(cv.fit)

library(penalized)
cvl.fit <- cvl(Y, penalized = X,model = "logistic",  lambda1 = grid)
#the fit on the full data (a penfit object)
cvl.fit$fullfit
#coefficients
cvl.fit$cvl

pen.fit <- penalized(Y, penalized = X,model = "logistic",  lambda1 = 1, steps = 100)



