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
  beta.hat.matrix = beta.hat.matrix[-1,]
  allzero.fun = function(beta){
    return(all(beta==0))
  }
  allzero.indicator = apply(beta.hat.matrix,2,allzero.fun)
  N = length(grid)
  w
  falses = which(!allzero.indicator )
  
  indicator = 1
  i = length(falses)
  while(indicator != 0 & i >= 1 ){
    j = falses[i]
    beta_j = beta.hat.matrix[,j]
    grid_j = grid[j]
    indicator.vec = sapply((j+1):N,function(k,x) infinity.norm(x[,k]-beta_j)/(grid[k]+grid_j)-C.bar<=0,
                           x = beta.hat.matrix)
    indicator = as.numeric(all(indicator.vec))
    i = i-1
    
  }  
  
  # j=N   
#   while(indicator != 0 & j > 1 ){
#     
#     beta_j = beta.hat.matrix[,j]
#     grid_j = grid[j]
#     indicator.vec = sapply(j:N,function(i,x) infinity.norm(x[,i]-beta_j)/(grid[i]+grid_j)-C.bar<=0,
#                            x = beta.hat.matrix)
#     indicator = as.numeric(all(indicator.vec))
#     j = j-1
#     
#   }
  
  lambda.hat = grid[j]
  return(lambda.hat)
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

fpfn.result = function(beta.star, coeffs, elaps, trial, end){
  fp.counts = fp.test(beta.star, coeffs)
  fn.counts = fn.test(beta.star, coeffs)
  end[trial,1] = fp.counts
  end[trial,2] = fn.counts
  end[trial,3] = fp.counts + fn.counts
  end[trial,4] = elaps  
  return(end)
}

nrow = 100; ncol_vec = c(100, 200, 500); 
sparsity_vec = c(6, 10, 15)
rho_vec = c(0, 0.25, 0.5)
C.bar=c(1.5)
for(sparsity in sparsity_vec){
  for(ncol in ncol_vec){
    for(rho in rho_vec){
      for(cbar in C.bar){
        beta.star = beta_star(ncol, sparsity)
        
        lambda.max=10*sqrt(log(ncol)/nrow)
        grid = seq(0.0001*lambda.max,lambda.max,length.out=500)
        simulations = 200
        end = matrix(0, simulations, 4)
        colnames(end) = c("fps", "fns", "jw.hd",  "jw.est.time") 
        
        end2 = matrix(0, simulations, 4)
        colnames(end2) = c("fps.cv", "fns.cv", "cv.hd", "cv.est.time")
        
        end3 = matrix(0, simulations, 4)
        colnames(end3) = c("fps.tr", "fns.tr", "hd.tr", "est.time.tr")
        
        for(trial in 1:simulations){
          try(
            {
              cat("n= ", nrow, "p= ", ncol, "rho= ", rho, "s= ", sparsity, "\n")
              cat("trial ", trial, " starts and Cbar = ", cbar, "\n")
              #define random design matrix
              set.seed(trial)
              
              desM = design_matrix(nrow, ncol, rho, beta.star)
              data = desM$data
              
              X.matrix = data[,-1]
              X.matrix = as.matrix(X.matrix)
              
              X.support = desM$suppMatrix
              X.support = as.matrix(X.support)
              
              #start time for alg
              t1=proc.time()
              lambda.hat = algorithm_new(data, grid, cbar)
              #end time for alg
              data=as.matrix(data)
              X=data[,-1]
              Y=data[,1]
              
              #find coefficients generated by using jw.lambda
              result = coeff.result(X, Y, lambda.hat)
              t2 = proc.time()
              
              ####   beta.hat will be thresholded by 3C.bar*lambda.hat
              beta.thresholding = beta_threshold(result[-1], factor=3, cbar, lambda.hat)
              t3 = proc.time()
              #store our elaps time 
              jw.elaps = t2-t1
              elaps.tr = t3-t1
              #lambda.hat
              
              #find fp counts, fn counts, fp+pn counts, estimate time by using jw.lambda and store them in end matrix
              end = fpfn.result(beta.star, result[-1], jw.elaps[3], trial, end)
              end3 = fpfn.result(beta.star, beta.thresholding, elaps.tr[3], trial, end3)
              
              
              #start time for cv
              t4=proc.time()
              cv.out = cv.glmnet(X,Y, family = "binomial", alpha = 1,lambda = grid)
              lambda.cv = cv.out$lambda.min
              
              
              #find coefficients generated by using cv.lambda
              result.cv = coeff.result(X, Y, lambda.cv)
              #stop time for cv 
              t5=proc.time()
              cv.elaps = t5-t4
              
              #find fp counts, fn counts, fp+pn counts, estimate time by using cv.lambda and stor them in end2 matrix
              end2 = fpfn.result(beta.star, result.cv[-1], cv.elaps[3], trial, end2)
              cat("trial ", trial, " ends and fpfn = ", end[trial, ], "\n")
              cat("trial ", trial, " ends and fpfn.cv = ", end2[trial, ], "\n")
              cat("trial ", trial, " ends and fpfn.tr = ", end3[trial, ], "\n")
            }
          )
        }
        
        
        end.matrix = data.frame(end, end2, end3)
        myfilename = paste0("r", nrow, "c", ncol, "o", rho, "s", sparsity, "Cbar" , cbar, ".csv")
        write.csv(end.matrix, file = myfilename)
      }
    }
  }
}
