
rm(list=ls())
library("mvtnorm")
library("glmnet")
library("methods")

expit=function(x){
  exp(x)/(1+exp(x))
}




# beta_star = function(p,sparsity,scale=c(-1,1)){
#   beta = rep(0,p)
#   index.zero = sample(1:p, p-sparsity, replace=FALSE)
#   beta[-index.zero] = sample(scale, sparsity, replace=TRUE, prob=c(0.5, 0.5))
#   return(beta)
# }

beta_star = function(p, sparsity, mu.beta, sigma = 0.02){
  beta = rep(0,p)
  index.zero = sample(1:p, p-sparsity, replace=FALSE)
  betas = log(rnorm(sparsity, mu.beta, 0.02))
  beta[-index.zero] = betas
  return(beta)
}

Ys = function(beta.star, X.support, support){
  #  dim(X.support) 
  logit.qs = NULL
  if(sum(support) == 1){
    nrow = length(X.support)
    ones = rep(1, nrow)
    u = ones %*% (X.support * beta.star[support]) / nrow
    logit.qs = u + X.support * beta.star[support]
  }else{
    nrow = nrow(X.support)
    ones = rep(1, nrow)
    u = ones %*% (X.support %*% beta.star[support]) / nrow   
    logit.qs = as.numeric(u) + X.support %*% beta.star[support]
  }
  
  
  #length(X.support); length(beta.star[support]); length(ones)
  
  qs = expit(logit.qs)
  ys = rbinom(nrow, 1, qs)
  return(ys)
}

design_matrix.A = function(nrow, ncol, beta.star){
  
  rho = 0
  support = (beta.star != 0)
  #length(beta.star)
  X = matrix(0, nrow, ncol) 
  #dim(X)
  
  Sigma=matrix(0,ncol,ncol)
  Sigma=(1-rho) * diag(ncol) + rho * (rep(1,ncol)) %o% (rep(1,ncol))
  X=rmvnorm(nrow,mean=rep(0,ncol),sigma=Sigma)
  #dim(X)
  X.support = X[, support]
  #dim(X.support)
  X.supp_comple=X[, !(support)]
  
  #generate ys 
  Y = Ys(beta.star, X.support, support)
  return(list( data = data.frame(Y,X), suppMatrix = X.support, supp_comMarix=X.supp_comple))
}

design_matrix.B = function(nrow, ncol, beta.star){
  support = (beta.star != 0)
  X = matrix(0, nrow, ncol) 
  Sigma=matrix(0,ncol,ncol)
  
  for(j in 1:ncol){
    Sigma[1:ncol %% 10 == j%%10, j] = 0.5
  }
  #   isSymmetric(Sigma)
  #   dim(X)
  #   dim(Sigma)
  X=rmvnorm(nrow,mean=rep(0,ncol),sigma=Sigma, method = "svd")
  X.support = X[,support]
  X.supp_comple=X[,!(support)]
  
  #generate ys 
  Y = Ys(beta.star,X.support, support)
  return(list( data = data.frame(Y,X), suppMatrix = X.support, supp_comMarix=X.supp_comple))
}


design_matrix.C = function(nrow, ncol, beta.star){
  support = (beta.star != 0)
  X = matrix(0, nrow, ncol) 
  Sigma=matrix(0,ncol,ncol)
  for(i in 1:ncol){
    for(j in 1:ncol){
      Sigma[i, j] = 0.9^(abs(i-j))
    }
  }
  #isSymmetric(Sigma)
  X=rmvnorm(nrow,mean=rep(0,ncol),sigma=Sigma)
  X.support = X[,support]
  X.supp_comple=X[,!(support)]
  
  #generate ys 
  Y = Ys(beta.star, X.support, support)
  return(list( data = data.frame(Y,X), suppMatrix = X.support, supp_comMarix=X.supp_comple))
}


design_matrix.D = function(nrow, ncol, beta.star){
  support = (beta.star != 0)
  X = matrix(0, nrow, ncol) 
  Sigma=matrix(0,ncol,ncol)
  for(i in 1:ncol){
    for(j in 1:ncol){
      Sigma[i, j] = 0.99^(abs(i-j))
    }
  }
  #isSymmetric(Sigma)
  X=rmvnorm(nrow,mean=rep(0,ncol),sigma=Sigma)
  X.support = X[,support]
  X.supp_comple=X[,!(support)]
  
  #generate ys 
  Y = Ys(beta.star, X.support, support)
  return(list( data = data.frame(Y,X), suppMatrix = X.support, supp_comMarix=X.supp_comple))
}

# 
# algorithm = function(data, grid, C.bar){
#   data=as.matrix(data)
#   X=data[,-1]
#   Y=data[,1]
#   glmnet.fit = glmnet(X, Y, family="binomial", alpha = 1, lambda = grid)
#   
#   beta.hat.matrix = coef(glmnet.fit, s = grid)
#   N = length(grid)
#   j = N
#   indicator = indicFunc(beta.hat.matrix, grid, C.bar, N, N)
#   while(indicator != 0 & j > 1){
#     j = j - 1
#     for (k in j:N){ 
#       indicator = indicFunc(beta.hat.matrix, grid, C.bar, j, k)
#       if(indicator == 0){
#         break
#       }
#     }
#   }
#   lambda.hat = grid[j]
#   return(lambda.hat)
# }



algorithm_new = function(data, grid, C.bar){
  data=as.matrix(data)
  X=data[,-1]
  Y=data[,1]
  glmnet.fit = glmnet(X, Y, family="binomial", alpha = 1, lambda = grid)
  beta.hat.matrix = as.matrix(coef(glmnet.fit, s = grid))
  beta.hat.matrix = data.frame(beta.hat.matrix)
  beta.hat.matrix = beta.hat.matrix[-1, ]
  allzero.fun = function(beta){
    return(all(beta==0))
  }
  allzero.indicator = apply(beta.hat.matrix, 2, allzero.fun)
  N = length(grid)
  falses = which(!allzero.indicator)
  # j=N
  indicator = 1
  i = length(falses)
  
  while(indicator != 0 & i >= 1){
    j = falses[i]
    beta_j = beta.hat.matrix[, j]
    grid_j = grid[j]
    indicator.vec = sapply((j + 1) : N,function(k, x) infinity.norm(x[, k]-beta_j)/(grid[k]+grid_j)-C.bar <= 0,
                           x = beta.hat.matrix)
    indicator = as.numeric(all(indicator.vec))
    i = i - 1
  }
  lambda.hat = grid[j]
  return(lambda.hat)
}


AIC_fun = function(data, grid){
  data=as.matrix(data)
  X=data[,-1]
  Y=data[,1]
  nrow = nrow(X)
  N = length(grid)
  glmnet.fit = glmnet(X, Y, family="binomial", alpha = 1, lambda = grid)
  beta.hat.matrix = as.matrix(coef(glmnet.fit, s = grid))
  beta.hat.matrix = data.frame(beta.hat.matrix)
  beta.hat.matrix = as.matrix(beta.hat.matrix)
  X_intercept = cbind(rep(1, nrow), X)
  AIC.values = sapply(1:N, function(j, x, y, coef){
    temp1 = as.vector(x %*% coef[, j])
    L_j = sum(y * temp1 - log(1 + exp(temp1)))
    nonzeros_j = sum(coef[-1, j] != 0)
    aic.vals = -2 * L_j + 2 * nonzeros_j
    return(aic.vals)
  }, x = X_intercept, y = Y, coef = beta.hat.matrix)
  index.min = which(AIC.values == min(AIC.values))
  index.lambda_hat = index.min[1]
  lambda.hat = grid[index.lambda_hat]
  return(lambda.hat)
}



BIC_fun = function(data, grid){
  data=as.matrix(data)
  X=data[,-1]
  Y=data[,1]
  nrow = nrow(X)
  N = length(grid)
  glmnet.fit = glmnet(X, Y, family="binomial", alpha = 1, lambda = grid)
  beta.hat.matrix = as.matrix(coef(glmnet.fit, s = grid))
  beta.hat.matrix = data.frame(beta.hat.matrix)
  beta.hat.matrix = as.matrix(beta.hat.matrix)
  X_intercept = cbind(rep(1, nrow), X)
  BIC.values = sapply(1:N, function(j, x, y, coef){
    temp1 = as.vector(x %*% coef[, j])
    L_j = sum(y * temp1 - log(1 + exp(temp1)))
    nonzeros_j = sum(coef[-1, j] != 0)
    bic.vals = -2 * L_j + log(nrow) * nonzeros_j
    return(bic.vals)
  }, x = X_intercept, y = Y, coef = beta.hat.matrix)
  index.min = which(BIC.values == min(BIC.values))
  index.lambda_hat = index.min[1]
  lambda.hat = grid[index.lambda_hat]
  return(lambda.hat)
}



infinity.norm=function(x){
  return(max(abs(x)))
}

# 
# indicFunc = function(beta.matrix, lambdas, C.bar, j, k){
#   beta.diff = beta.matrix[-1, j] - beta.matrix[-1,k]
#   lambda.sum = lambdas[j] + lambdas[k]
#   func = newnorm(beta.diff)/lambda.sum - C.bar
#   indicator = 0
#   if(func <= 0){
#     indicator = 1
#   }
#   return (indicator)
# }
# 
# 
# l_infty=function(beta.star,beta.hat){
#   beta.diff=beta.hat-beta.star
#   return(newnorm(beta.diff))
# }

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

power = function(base, test){
  index.nonzeros = which(base != 0)
  length.nonzeros = length(index.nonzeros)
  true.discover = ifelse(test[index.nonzeros] != 0, 1, 0)
  return(sum(true.discover) / length.nonzeros)
}


oneminus.fdr = function(base, test){
  index.zeros = which(base == 0)
  length.nonzeros.test = length(which(test != 0))
  non.zeros = ifelse(test[index.zeros] != 0, 1, 0)
  if(length.nonzeros.test != 0){
    oneminus_fdr = 1 - (sum(non.zeros) / length.nonzeros.test)
  }else{
    oneminus_fdr = 1
  }
  return(oneminus_fdr)
}

fpfn.result = function(beta.star, coeffs, elaps, trial, end){
  fp.counts = fp.test(beta.star, coeffs)
  fn.counts = fn.test(beta.star, coeffs)
  power.value = power(beta.star, coeffs)
  oneminus_fdr = oneminus.fdr(beta.star, coeffs)
  end[trial,1] = power.value
  end[trial,2] = oneminus_fdr
  end[trial,3] = power.value + oneminus_fdr
  end[trial,4] = fp.counts + fn.counts  
  end[trial, 5] = elaps
  return(end)
}


mu.beta.lows = c(log(1.35), log(1.75)) #low dimension corresponding to n = 1000
nrow = 200 # n = 200 is high-dimension, n = 1000 is low-dimension

ncol = 500;

sparsity_vec = c(5, 10, 20)
#sparsity = 5

cbar = 1.5
simulations = 200

k = 1



mu.beta = log(2.5)
sparsity = 5
trial = 1
for(mu.beta in mu.beta.lows){
  for(sparsity in sparsity_vec){  
    
    set.seed(k)
    beta.star = beta_star(ncol, sparsity, mu.beta)
    
    lambda.max=20*sqrt(log(ncol)/nrow)
    grid = seq(0.0001*lambda.max,lambda.max,length.out=500)
    
    end.av = matrix(0, simulations, 5)
    colnames(end.av) = c("oneminus.av", "power.av", "hd.av",  "hd2.av", "elasp")  #hd2 = fp.counts + fn.counts
    
    end.tr = matrix(0, simulations, 5)
    colnames(end.tr) = c("oneminus.tr", "power.tr", "hd.tr", "hd2.tr", "elasp")
    
    end.cv = matrix(0, simulations, 5)
    colnames(end.cv) = c("oneminus.cv", "power.cv", "hd.cv", "hd2.cv", "elasp")
    
    end.aic = matrix(0, simulations, 5)
    colnames(end.aic) = c("oneminus.aic", "power.aic", "hd.aic", "hd2.aic", "elasp")
    
    end.bic = matrix(0, simulations, 5)
    colnames(end.bic) = c("oneminus.bic", "power.bic", "hd.bic", "hd2.bic", "elasp")
    
    for(trial in 1:simulations){ 
      try( 
        { 
          #trial = 1
          #define random design matrix
          set.seed(trial + 100)
          desM = design_matrix.C(nrow, ncol, beta.star)
          
          cat("n= ", nrow, "p= ", ncol, "mu.beta= ", mu.beta, "s= ", sparsity, "desM.C and HD",  "\n") #HD is high dimension
          
          data = desM$data
          
          X.matrix = data[,-1]
          X.matrix = as.matrix(X.matrix)
          
          X.support = desM$suppMatrix
          X.support = as.matrix(X.support)
          
          #start time for alg
          t_start.av=proc.time()
          lambda.hat = algorithm_new(data, grid, cbar)
          #end time for alg
          data=as.matrix(data)
          X=data[,-1]
          Y=data[,1]
          
          #find coefficients generated by using jw.lambda
          result = coeff.result(X, Y, lambda.hat)
          t_end.av = proc.time()
          
          ####   beta.hat will be thresholded by 3C.bar*lambda.hat
          beta.thresholding = beta_threshold(result[-1], factor=3, cbar, lambda.hat)
          t_end.tr = proc.time()
          #store our elaps time 
          elaps.av = t_end.av - t_start.av
          elaps.tr = t_end.tr - t_start.av
          #lambda.hat
          
          #find oneminus fdr, average power, hd1, fp+pn counts, estimate time by using jw.lambda and store them in end matrix
          end.av = fpfn.result(beta.star, result[-1], elaps.av[3], trial, end.av)
          
          end.tr = fpfn.result(beta.star, beta.thresholding, elaps.tr[3], trial, end.tr)
          
          #start time for cv
          t_start.cv = proc.time()
          cv.out = cv.glmnet(X,Y, family = "binomial", alpha = 1,lambda = grid)
          lambda.cv = cv.out$lambda.min
          
          
          #find coefficients generated by using cv.lambda
          result.cv = coeff.result(X, Y, lambda.cv)
          #stop time for cv 
          t_end.cv = proc.time()
          elaps.cv = t_end.cv - t_start.cv
          
          #find oneminus fdr, average power, hd1, fp+pn counts, estimate time by using cv.lambda and stor them in end2 matrix
          end.cv = fpfn.result(beta.star, result.cv[-1], elaps.cv[3], trial, end.cv)
          ### AIC methods choosing optimal tuning parameter   
          t_start.aic = proc.time()
          lambda.aic = AIC_fun(data, grid)
          result.aic = coeff.result(X, Y, lambda.aic)
          t_end.aic = proc.time()
          elaps.aic = t_end.aic - t_start.aic
          end.aic = fpfn.result(beta.star, result.aic[-1], elaps.aic[3], trial, end.aic)
          ### BIC methods
          t_start.bic = proc.time()
          lambda.bic = BIC_fun(data, grid)
          result.bic = coeff.result(X, Y, lambda.bic)
          t_end.bic = proc.time()
          elaps.bic = t_end.bic - t_start.bic
          end.bic = fpfn.result(beta.star, result.bic[-1], elaps.bic[3], trial, end.bic)
          
          cat("designM(Cl",  ") trial ", trial, " ends. ",  "(av)h1 and h2 are ", end.av[trial, 3]," and ", end.av[trial, 4], "\n")
          cat("designM(Cl",  ") trial ", trial, " ends. ",  "(tr)h1 and h2 are ", end.tr[trial, 3], " and ", end.tr[trial, 4], "\n")
          cat("designM(Cl",  ") trial ", trial, " ends. ",  "(cv)h1 and h2 are = ", end.cv[trial, 3], " and ", end.cv[trial, 4], "\n")
          cat("designM(Cl",  ") trial ", trial, " ends. ",  "(aic)h1 and h2 are = ", end.aic[trial, 3], " and ", end.aic[trial, 4], "\n")
          cat("designM(Cl",  ") trial ", trial, " ends. ",  "(bic)h1 and h2 are = ", end.bic[trial, 3], " and ", end.bic[trial, 4], "\n")
          
        }
      )
    }
    
    
    end.matrix = data.frame(end.av, end.tr, end.cv, end.aic, end.bic)
    myfilename = paste0("Cl", "r", nrow, "c", ncol,  "u", mu.beta, "s", sparsity, ".csv")
    write.csv(end.matrix, file = myfilename)
    
    
    
    k = k+1
    
    
  }    
  
  
}


