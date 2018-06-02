set.seed(1)
#install.packages("mvtnorm")
library("mvtnorm")
library(glmnet)


expit=function(x){
  exp(x)/(1+exp(x))
}


nrow = 50; ncol = 100; beta = runif(100)
threshold = 0.9
rho = 1
beta.star = beta_star(beta, threshold)
desM = design_matrix(nrow, ncol, rho, beta.star)
names(desM)

factor = 100
data = desM$data
X.matrix = data[,-1]
X.matrix = as.matrix(X.matrix)

X.support = desM$suppMatrix
dim(X.support)
X.support = as.matrix(X.support)

W.star = diag_new(X.matrix, beta.star)
X.supp_comple=desM$supp_comMarix ###

##here need help debugging one_Minus_b() to deal with solve() error
lambda.max = findLambdaMax(X.matrix, X.support, W.star, X.supp_comple)
grid = gridFunc(lambda.max, factor)
C.bar = findCbar(X.support, W.star)

lambda.hat = algorithm(data, grid, C.bar)
##design_matrix cannot return X.support_comple matrix




design_matrix = function(nrow, ncol, rho, beta.star){
  support = (beta.star != 0)
  #  support = beta.star[support.index]
  X = matrix(0, nrow, ncol) ##initialaztion
  Sigma=matrix(0,ncol,ncol)
  Sigma=rho*diag(ncol)+(1-rho)*(rep(1,ncol))%o%(rep(1,ncol))
  X=rmvnorm(nrow,mean=rep(0,ncol),sigma=Sigma)
  X.support = X[,support]
  X.supp_comple=X[,-(support)]
  temp=runif(nrow)
  Y = as.numeric(temp<=expit(X %*% beta.star))
  return(list( data = data.frame(Y,X),suppMatrix = X.support, supp_comMarix=X.supp_comple))
  #why cannot return supp_comMatix
}



algorithm = function(data, grid, C.bar){
  data=as.matrix(data)
  X=data[,-1]
  Y=data[,1]
  glmnet.fit = glmnet(X, Y, family="binomial", alpha = 1, lambda = grid)
  
  beta.hat.matrix = coef(glmnet.fit, s = grid)
  #######here I need to transform beta.hat.matrix to numeric matrix
  
  str(beta.hat.matrix)
  dim(beta.hat.matrix)
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




gridFunc = function(lambda.max, factor){
  return(c(1:factor)/factor * lambda.max)
}
findCbar = function(X.support, W.star){
  a.value = a(X.support, W.star)
  C.min = Cmin(X.support, W.star)
  return( Cbar(a.value, C.min)) 
}
findLambdaMax = function(X.matrix, X.support, W.star, X.supp_comple){
  C.min = Cmin(X.support, W.star) + 1## choose bigger than c.min
  D.max = Dmax(X.matrix)
  b.value = 1 - one_Minus_b(X.support, W.star, X.supp_comple) + 19 #choose
  L.s = L(X.support)
  return(lambda_max(b.value, C.min, D.max, L.s ))
}

support_matrix = function(design.matrix, support){
  support.Matrix = design.matrix[, support]
}
beta_star = function(beta, threshold){
  beta[which(beta<threshold)] = 0
  return(beta)
}
##Cbar greater or equal to Cbar()
Cbar = function(a.value, C.min){
  return(3*a.value / (2 * C.min))
}
Cmin = function(X.support, W.star){
  C = t(X.support) %*% W.star %*% X.support/nrow(W.star)
  return(min(eigen(C)$values)) 
}
Dmax = function(X.matrix){
  D = t(X.matrix) %*% X.matrix / nrow(X.matrix)
  return(max(eigen(D)$values))
}
#a should be at least greater the value returned by a() function
a = function(X.support, W.star){
  matrixA = square_matrix(X.support, W.star)
  matrixA.inv = solve(t(matrixA)) ################need debug
  numerator = norm(matrixA.inv, "I")
  denorminator = norm(matrixA.inv, "E")
  return(numerator/denorminator)
}

#1-b should be at least greater than the one_Minus_b() returned value
one_Minus_b = function(X.support, W.star, X.supp_comple){
  matrixB = square_matrix(X.support, W.star)
  matrixB.inv = solve(t(matrixB)) ###############need debug
  matrixBb = matrixB.inv %*% t(X.support) %*% W.star %*% X.supp_comple
  return(norm(matrixBb, "I"))
}

#L should be at least greater than the returned L() value+
L = function(X.matrix){
  return (max(X.matrix))
}

lambda_max = function(b, C.min, D.max, Ls ) {
  return (b*C.min^2)/(100*(2-b)*D.max*Ls)
}

newnorm=function(x){
  return(max(abs(x)))
}


indicFunc = function(beta.matrix, lambdas, C.bar, j, k){
  beta.diff = beta.matrix[, j] - beta.matrix[,k]
  lambda.sum = lambdas[j] + lambdas[k]
  func = newnorm(beta.diff)/lambda.sum - C.bar
  indicator = 0
  if(func <= 0){
    indicator = 1
  }
  return (indicator)
}

diag_new = function(matrixX, beta){
  Wstar = diag(rep(1, nrow(matrixX)))
  for (i in 1:nrow) {
    Wstar[i,i] = Wfunc(matrixX[i,], beta)
  }
  return(Wstar)
}

Wfunc = function(row, beta){
  return(exp(row %*% beta)/(1+exp(row %*% beta))^2)
}

square_matrix = function(matrixX, transposer){
  Smatrix = t(matrixX) %*% transposer %*% matrixX
  return(Smatrix)
}