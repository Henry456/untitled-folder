rm(list=ls())
library("mvtnorm")
library("glmnet")
library("methods")

expit=function(x){
  exp(x)/(1+exp(x))
}

design_matrix = function(nrow, ncol, rho, beta.star){
  support = (beta.star != 0)
  X = matrix(0, nrow, ncol) 
  Sigma=matrix(0,ncol,ncol)
  Sigma=(1-rho) * diag(ncol) + rho * (rep(1,ncol)) %o% (rep(1,ncol))
  X=rmvnorm(nrow,mean=rep(0,ncol),sigma=Sigma)
  X.support = X[,support]
  X.supp_comple=X[,!(support)]
  temp=runif(nrow)
  Y = as.numeric(temp<=expit(X %*% beta.star))
  return(list( data = data.frame(Y,X), suppMatrix = X.support, supp_comMarix=X.supp_comple))
}

beta_star = function(p,sparsity,scale=c(-1,1)){
  beta=rep(0,p)
  index.zero=sample(1:p,p-sparsity,replace=FALSE)
  beta[-index.zero]=sample(scale,sparsity,replace=TRUE,prob=c(0.5,0.5))
  return(beta)
}

infinity.norm=function(x){
  return(max(abs(x)))
}

fp.test = function(base, test){
  index.zeros = which(base == 0)
  non.zeros = ifelse(test[index.zeros] != 0, 1, 0) #only contains 1s, 0s.
  return (sum(non.zeros))
}

fn.test = function(base, test){
  index.nonzeros = which(base != 0)
  zeros = ifelse(test[index.nonzeros] == 0, 1, 0)
  return (sum(zeros))
}

coeff.result = function(X, Y, lambda.hat){
  glmnet.fit = glmnet(X, Y, family="binomial", alpha = 1, lambda = lambda.hat)
  result=coef(glmnet.fit)
  result = as.vector(result)
  return(result)
}

beta_threshold = function(beta,factor=3,C.bar,lambda){
  beta[abs(beta)<factor*C.bar*lambda]=0
  return(beta)
}

fpfn.result = function(beta.star, coeffs, elaps, trial, end){
  fp.counts = fp.test(beta.star, coeffs)
  fn.counts = fn.test(beta.star, coeffs)
  end[trial,1] = fp.counts
  end[trial,2] = fn.counts
  end[trial,3] = fp.counts + fn.counts
  end[trial,4] = elaps  
  return(end)
}


nrow = 50; ncol = 50
sparsity = 4
rho = 0
C.bar=c(1.5)


beta.star = beta_star(ncol, sparsity)

lambda.max=10*sqrt(log(ncol)/nrow)
grid = seq(0.0001*lambda.max,lambda.max,length.out=500)
simulations = 1
end = matrix(0, simulations, 4)
colnames(end) = c("fps", "fns", "jw.hd",  "jw.est.time") 

end2 = matrix(0, simulations, 4)
colnames(end2) = c("fps.cv", "fns.cv", "cv.hd", "cv.est.time")

end3 = matrix(0, simulations, 4)
colnames(end3) = c("fps.tr", "fns.tr", "hd.tr", "est.time.tr")
set.seed(trial)

desM = design_matrix(nrow, ncol, rho, beta.star)
data1 = desM$data

allzero.fun = function(beta){
  return(all(beta==0))
}


#algorithm_new = function(data1, grid, C.bar){
  
  data=as.matrix(data1)
  X=data[,-1]
  Y=data[,1]
  glmnet.fit = glmnet(X, Y, family="binomial", alpha = 1, lambda = grid)
  beta.hat.matrix = as.matrix(coef(glmnet.fit, s = grid))
  beta.hat.matrix = data.frame(beta.hat.matrix)
  beta.hat.matrix = beta.hat.matrix[-1,]

  
  allzero.indicator = apply(beta.hat.matrix,2,allzero.fun)
  N = length(grid)
  
  falses = which(!allzero.indicator )
  beta.hat.matrix.falses = beta.hat.matrix[, falses]
  
  indicators = sapply ( length(falses):1, 
                     function(i, x) {
                        indicator.vec = sapply(falses[i]:N,

                                           function(k,x) {
                                 
                                             indicator.vect = infinity.norm(x[,k]-beta.hat.matrix.falses[,i])/(grid[k]+grid[falses[i]])-C.bar<=0
                                             
                                             
                                           }, 
                        
                                           x = beta.hat.matrix
                                    )
                        indicator = as.numeric(all(indicator.vec))
                        return(indicator)
                      },
                      
                      x = beta.hat.matrix.falses 
                ) 
  i = which(indicators == F)[1] - 1 #The first False - 1 = the last True 
  j = falses[i]
  lambda.hat = grid[j]
#  return(lambda.hat)
#


