library(glmnet)
#low dimension sparse linear regression
library(ISLR)
#fix(Hitters)
Hitters=na.omit(Hitters)
#use the low dimentional data called Hitters

#construct a matrix x which contains the features
#and y which we want to predict
x <- model.matrix(Salary~., Hitters)[, -1]
#fix(x)
y <- Hitters$Salary
#grab some data from the sample above as training data
#and the leftover as the test data
train=sample(1:nrow(x), nrow(x)/2)
test=(-train)
y.test=y[test]

#to construct a lasso model with 100 lambdas made manually
#which will be input in the glmnet
#later we will use cross-validation method to choose the best lambda
grid <- 10^seq(10,-2,length=100)
lasso.mod <- glmnet(x[train,], y[train], alpha=1, lambda=grid)

#cross-validation to choose the best lambda, which contains
#the lambda.min and lambda.1se, based on the trainning set
set.seed(1)
cv.out <- cv.glmnet(x[train,], y[train], alpha=1)
#plot(cv.out)
bestlam <- cv.out$lambda.min
#based the previous lass model to predict over test data set
lasso.pred <- predict(lasso.mod, s=bestlam, newx=x[test,])
#calculate the mean square error
mean((lasso.pred-y.test)^2)

out <- glmnet(x, y, alph = 1, lambda = grid)
lasso.coef=predict(out, type="coefficients", s=bestlam)[1:20,]
#to achieve the coefficients of all the features
lasso.coef
#to achieve non-zero coefficients of features
lasso.coef[lasso.coef!=0]
#eq, coef(out, s=bestlam)


#(2)high-dimensional sparse linear regression
#clean the data
riboflavin <- read.csv("riboflavin.csv")
cnames <- riboflavin[, 1]
#fix(riboflavin)
#Need to convert to matrix by using data.matrix (important)
riboflavin <- data.matrix(riboflavin)
#fix(riboflavin)

riboflavin <- t(riboflavin)
#fix(riboflavin)
riboflavin <- riboflavin[-1, ]
#fix(riboflavin)
colnames(riboflavin) <- cnames
#fix(riboflavin)
riboflavin <- na.omit(riboflavin)


x <- riboflavin[, -1]
y <- riboflavin[, 1]


fit <- glmnet(x = x, y = y, alpha = 1)
plot(fit)
set.seed(1)
fit.cv <- cv.glmnet(x = x, y = y, alpha = 1)


plot(fit.cv)
bestlam <- fit.cv$lambda.min
lasso.coef=predict(fit, type="coefficients", s=bestlam)[1:4089,]
#coef(fit, s=bestlam)
lasso.coef[lasso.coef!=0]
#> length(lasso.coef[lasso.coef!=0])
#[1] 41
#which indicate that the number
#of non-zero  coefficients are 41

