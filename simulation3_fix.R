library("mvtnorm")
library(glmnet)

design_matrix = function(nrow, ncol, rho, beta.star){
  support = (beta.star != 0)
  X = matrix(0, nrow, ncol) 
  Sigma=matrix(0,ncol,ncol)
  Sigma=rho*diag(ncol)+(1-rho)*(rep(1,ncol))%o%(rep(1,ncol))
  X=rmvnorm(nrow,mean=rep(0,ncol),sigma=Sigma)
  X.support = X[,support]
  X.supp_comple=X[,!(support)]
  temp=runif(nrow)
  Y = as.numeric(temp<=expit(X %*% beta.star))
  return(list( data = data.frame(Y,X), suppMatrix = X.support, supp_comMarix=X.supp_comple))
}


algorithm = function(data, grid, C.bar){
  data=as.matrix(data)
  X=data[,-1]
  Y=data[,1]
  glmnet.fit = glmnet(X, Y, family="binomial", alpha = 1, lambda = grid)
  
  beta.hat.matrix = coef(glmnet.fit, s = grid)
  N = length(grid)
  j = N
  indicator = indicFunc(beta.hat.matrix, grid, C.bar, N, N)
  while(indicator != 0 & j > 1){
    j = j - 1
    for (k in j:N){ 
      indicator = indicFunc(beta.hat.matrix, grid, C.bar, j, k)
      if(indicator == 0){
        break
      }
    }
  }
  lambda.hat = grid[j]
  return(lambda.hat)
}


beta_star = function(p,sparsity,scale=c(-1,1)){
  beta=rep(0,p)
  index.zero=sort(sample(1:p,p-sparsity,replace=FALSE))
  beta[-index.zero]=sample(scale,sparsity,replace=TRUE,prob=c(0.5,0.5))
  return(beta)
}

newnorm=function(x){
  return(max(abs(x)))
}


indicFunc = function(beta.matrix, lambdas, C.bar, j, k){
  beta.diff = beta.matrix[-1, j] - beta.matrix[-1,k]
  lambda.sum = lambdas[j] + lambdas[k]
  func = newnorm(beta.diff)/lambda.sum - C.bar
  indicator = 0
  if(func <= 0){
    indicator = 1
  }
  return (indicator)
}


l_infty=function(beta.star,beta.hat){
  beta.diff=beta.hat-beta.star
  return(newnorm(beta.diff))
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

nrow = 300; ncol = 100; 
sparsity=6
rho = 1
beta.star = beta_star(ncol, sparsity)

lambda.max=10*sqrt(log(ncol)/nrow)
grid = seq(0.0001*lambda.max,lambda.max,length.out=500)

C.bar=c(6,3,1.5,2.6,2,1.8,1.2,1,0.9,1.3)

simulations = 2000
end = matrix(0, simulations, 4)
names(end) = c("fps", "fn", "est.time", "cv.compare.time")
for(trial in 1:simulations){
  
  #define random design matrix
  set.seed(trial)
  
  desM = design_matrix(nrow, ncol, rho, beta.star)
  data = desM$data
  
  X.matrix = data[,-1]
  X.matrix = as.matrix(X.matrix)
  
  X.support = desM$suppMatrix
  X.support = as.matrix(X.support)
  
  #start time for alg
  t1=Sys.time()
  lambda.hat = algorithm(data, grid, C.bar[8])
  #end time for alg
  t2=Sys.time()
  
  #store our elaps time 
  jw.elaps = t2-t1
  end[trial,3] = jw.elaps
  
  #lambda.hat
  data=as.matrix(data)
  X=data[,-1]
  Y=data[,1]
  glmnet.fit = glmnet(X, Y, family="binomial", alpha = 1, lambda = lambda.hat)
  result=coef(glmnet.fit)
  result = as.vector(result)
  l.inftyerror=l_infty(beta.star,result[-1])
  l.inftyerror
  
  #find fp counts, fn counts and stor them in end matrix
  fp.counts = fp.test(beta.star, result[-1])
  fn.counts = fn.test(beta.star, result[-1])
  end[trial,1] = fp.counts
  end[trial,2] = fp.counts
  
  which(result[-1]!=0)
  which(beta.star!=0)
  
  t1=Sys.time()
  cv.out = cv.glmnet(X,Y, alpha = 1, lambda = grid)
  t2=Sys.time()
  cv.elaps = t2-t1 ####very small. Is this what you want me to test?
  end[trial, 4] = cv.elaps
}

head(end)
